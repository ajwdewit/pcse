# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
"""Routines for retrieving data from the PCSE demo database.
 
Implements the following functions:
    - fetch_cropdata()
    - fetch_sitedata()
    - fetch_soildata()
    - fetch_timerdata()

"""
import sys, os
import datetime as dt
import logging

import yaml

from ..exceptions import PCSEError
from ..base import WeatherDataContainer, WeatherDataProvider
from ..util import wind10to2, reference_ET, safe_float, check_date
from .. import settings


def fetch_crop_name(DBconn, crop):
    # Get crop name from crop table
    cursor = DBconn.cursor()
    cursor.execute("select crop_name from crop where crop_no=?", (crop,))
    row = cursor.fetchone()
    if row:
        return row.crop_name
    else:
        raise PCSEError(f"No crop_name found for crop_no={crop}")


def fetch_cropdata(DBconn, grid, year, crop):
    """Retrieve crop parameter values for given grid, year, crop from DB.
    
    Parameter values are pulled from tables 'crop_parameter_value'
    and 'variety_parameter_value'. Metadata is an SQLAlchemy metadata object.
    
    Returns a dictionary with WOFOST crop parameter name/value pairs.
    
    Note that the parameter names are defined in the function itself in order
    to distinguish scalar and table parameters. These definitions will need to
    be extended when additional parameters are to be retrieved.
    """
    
    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger(__name__)

    # Create initial dictionary 
    cropdata = {}
    cropdata["CRPNAM"] = fetch_crop_name(DBconn, crop)

    # Get crop variety from crop_calendar;
    cursor = DBconn.cursor()
    cursor.execute("select variety_no from crop_calendar where crop_no=? and grid_no=? and year=?", (crop, grid, year))
    rows = cursor.fetchall()
    variety = rows[0].variety_no

    # Define crop parameter values
    parameter_codes_sngl = ("CFET", "CVL", "CVO", "CVR", "CVS", "DEPNR", "DLC", 
                            "DLO", "DVSEND", "EFF", "IAIRDU", "IDSL", "KDIF", 
                            "LAIEM", "PERDL", "Q10", "RDI", "RDMCR", "RGRLAI", 
                            "RML", "RMO", "RMR", "RMS", "RRI", "SPA", "SPAN", "SSA", 
                            "TBASE", "TBASEM", "TDWI", "TEFFMX", "TSUM1", "TSUM2", 
                            "TSUMEM", "IOX")
    parameter_codes_mltp = ("AMAXTB", "DTSMTB", "FLTB", "FOTB", "FRTB", "FSTB", 
                            "RDRRTB", "RDRSTB", "RFSETB", "SLATB", "TMNFTB", 
                            "TMPFTB")
    
    # Pull single value parameters from CROP_PARAMETER_VALUE first
    sql = "select * from crop_parameter_value where crop_no=? and parameter_code=?"
    for paramcode in parameter_codes_sngl:
        r = cursor.execute(sql, (crop, paramcode))
        rows = r.fetchall()
        cropdata[paramcode] = float(rows[0].parameter_xvalue)

    # Pull array parameter values from CROP_PARAMETER_VALUE
    # note the change in the mask value and the use of "LIKE" in the SQL query
    sql = "select * from crop_parameter_value where crop_no=? and parameter_code like ?"
    for paramcode in parameter_codes_mltp:
        pattern = paramcode + r'%'
        r = cursor.execute(sql, (crop, pattern))
        values = []
        for row in r.fetchall():
            values.append(float(row.parameter_xvalue))
            values.append(float(row.parameter_yvalue))
        cropdata[paramcode] = values

    # Pull same parameter values from VARIETY_PARAMETER_VALUES
    # if they are defined for that variety.
    # Pull single value parameters first
    sql = "select * from variety_parameter_value where crop_no=? and variety_no=? and parameter_code=?"
    for paramcode in parameter_codes_sngl:
        r = cursor.execute(sql, (crop, variety, paramcode))
        rows = r.fetchall()
        if rows:
            cropdata[paramcode] = float(rows[0].parameter_xvalue)

    # pull array value parameters - note the change in the mask value and
    # the use of "LIKE" in the SQL query
    sql = "select * from variety_parameter_value where crop_no=? and variety_no=? and parameter_code like ?"
    for paramcode in parameter_codes_mltp:
        pattern = paramcode + r'%'
        r = cursor.execute(sql, (crop, variety, pattern))
        rows = r.fetchall()
        if rows:
            values = []
            for row in rows:
                values.append(float(row.parameter_xvalue))
                values.append(float(row.parameter_yvalue))
            cropdata[paramcode] = values

    cursor.close()

    # Make some specific changes for PCSE wofost 7.2 with regard to variables
    # SSA, KDIF and EFF. This is needed because the 7.2 code expects a
    # parameter array, while these parameters have been defined in CGMS as
    # single values. DVSI does not exist in CGMS, therefore set DVSI to zero.
    
    # SSA convert to SSATB:
    SSA = cropdata["SSA"]
    cropdata.update({"SSATB": [0, SSA, 2.0, SSA]})
    # KDIF convert to KDIFTB:
    KDIF = cropdata["KDIF"]
    cropdata.update({"KDIFTB": [0., KDIF, 2.0, KDIF]})
    # EFF convert to EFFTB
    EFF = cropdata["EFF"]
    cropdata.update({"EFFTB": [0., EFF, 40.0, EFF]})
    # DVSI set to 0
    cropdata.update({"DVSI":0})
    
    logger.info("Succesfully retrieved crop parameter values from database")
    return cropdata


def fetch_soildata(DBconn, grid):
    """Retrieve soil parameters for given grid from DB for a 1-layer soil.
    
    Retrieves soil_type_no from the table SOIL_TYPE and associated soil layers
    and soil physical data from tables SOIL_LAYERS and SOIL_PHYSICAL_GROUP. 
    
    Returns a dictionary with WOFOST soil parameter name/value pairs.
    """
    
    cursor = DBconn.cursor()

    soildata = {}
    # Select soil from the table SOIL_TYPE
    sql = "select * from soil_type where grid_no=?"
    r = cursor.execute(sql, (grid,))
    row = r.fetchone()
    soil_type_no = row.soil_type_no
    
    # Derive layers for given soil_type_no. This should return only one
    # layer, otherwise raise an error.
    sql = "select * from soil_layers where soil_type_no=? order by layer_no"
    r  = cursor.execute(sql, (soil_type_no,))
    rows = r.fetchall()
    if len(rows) == 0:
        msg = "No record found."
        raise PCSEError(grid, msg)
    elif len(rows) > 1:
        msg = ("Number of soil layers > 1. Not possible for unlayered " +
               "waterbalance module. Use 'fetch_soiltype_multilayer'") 
        raise PCSEError(grid, msg)
    else:
        soildata["RDMSOL"] = float(rows[0].thickness)
        soil_group_no = rows[0].soil_group_no
    
    # Retrieve soil physical properties for given layer for given soil
    # parameter codes: (wofost_parname, database_name)
    soil_parameters = [("CRAIRC", "CRITICAL_AIR_CONTENT"),
                       ("K0", "HYDR_CONDUCT_SATUR"),
                       ("SOPE", "MAX_PERCOL_ROOT_ZONE"),
                       ("KSUB", "MAX_PERCOL_SUBSOIL"),
                       ("SMFCF", "SOIL_MOISTURE_CONTENT_FC"),
                       ("SM0", "SOIL_MOISTURE_CONTENT_SAT"),
                       ("SMW", "SOIL_MOISTURE_CONTENT_WP")]
    # table_soil_pg = Table('soil_physical_group',metadata, autoload=True)
    sql = "select * from soil_physical_group where soil_group_no=? and parameter_code=?"
    for (wofost_soil_par, db_soil_par) in soil_parameters:
        r = cursor.execute(sql, (soil_group_no, db_soil_par))
        row = r.fetchone()
        if row is None:
            msg = "Parameter %s not found" % db_soil_par
            raise PCSEError(grid, msg)
        soildata[wofost_soil_par] = float(row.parameter_xvalue)

    return soildata


class AgroManagementDataProvider(list):
    """Class for providing agromanagement data from the CROP_CALENDAR table in a PCSE database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the grid_no column in the table
    :param crop_no: Integer id of crop, maps to the crop_no column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
           The campaign year refers to the year of the crop start. Thus for crops
           crossing calendar years, the start_date can be in the previous year as the
           harvest.
    
    Note that this AgroManagementDataProvider is only used for the internal PCSE database
    and not to be used for CGMS databases.
    """
    agro_management_template = """
          - {campaign_start_date}:
                CropCalendar:
                    crop_name: '{crop_name}'
                    variety_name: '{variety_name}'
                    crop_start_date: {crop_start_date}
                    crop_start_type: {crop_start_type}
                    crop_end_date: {crop_end_date}
                    crop_end_type: {crop_end_type}
                    max_duration: {duration}
                TimedEvents: null
                StateEvents: null
        """

    def __init__(self, DBconn, grid_no, crop_no, campaign_year):
        list.__init__(self)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.crop_name = fetch_crop_name(DBconn, self.crop_no)

        cursor = DBconn.cursor()
        sql = "select * from crop_calendar where grid_no=? and crop_no=? and year=?"
        r = cursor.execute(sql, (self.grid_no, self.crop_no, self.campaign_year))
        row = r.fetchone()
        if not row:
            msg = f"Failed deriving crop calendar for grid_no {self.grid_no}, crop_no {self.crop_no}, year {self.campaign_year}"
            raise PCSEError(msg)

        # Determine the start date.
        self.crop_start_date = check_date(row.crop_start_date)
        self.campaign_start_date = row.start_date

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        self.crop_start_type = str(row.crop_start_type).strip()
        if self.crop_start_type not in ["sowing","emergence"]:
            msg = "Unrecognized crop start type: %s" % self.crop_start_type
            raise PCSEError(msg)

        # Determine maximum duration of the crop
        self.max_duration = int(row.max_duration)

        # Determine crop end date/type and the end of the campaign
        self.crop_end_type = str(row.crop_end_type).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % row.end_type)
            raise PCSEError(msg)

        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
        else:
            self.crop_end_date = check_date(row.crop_end_date)

        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def _build_yaml_agromanagement(self):
        """Builds the YAML agromanagent string"""

        # We do not get a variety_name from the CGMS database, so we make one
        # as <crop_name>_<grid>_<year>
        variety_name = "%s_%s_%s" % (self.crop_name, self.grid_no, self.campaign_year)
        input = self.agro_management_template.format(campaign_start_date=self.campaign_start_date,
                                                     crop_name=self.crop_name,
                                                     variety_name=variety_name,
                                                     crop_start_date=self.crop_start_date,
                                                     crop_start_type=self.crop_start_type,
                                                     crop_end_date=self.crop_end_date,
                                                     crop_end_type=self.crop_end_type,
                                                     duration=self.max_duration
                                                     )
        return input

    def _parse_yaml(self, input):
        """Parses the input YAML string and assigns to self"""
        try:
            items = yaml.safe_load(input)
        except yaml.YAMLError as e:
            msg = "Failed parsing agromanagement string %s: %s" % (input, e)
            raise PCSEError(msg)
        self.extend(items)


def fetch_sitedata(DBconn, grid, year):
    """Retrieve site data from DB for given grid, year.
    
    Pulls sitedata from the PCSE database 'SITE' table,
    
    Returns a dictionary with site parameter name/value pairs.
    """

    cursor = DBconn.cursor()
    r = cursor.execute("select * from site where grid_no=? and year=?", (grid, year))
    row = r.fetchone()
    if row is not None:
        sitedata = {}
        sitedata['IFUNRN'] = float(row.ifunrn)
        sitedata['SSMAX'] = float(row.max_surface_storage)
        sitedata['NOTINF'] = float(row.not_infiltrating_fraction)
        sitedata['SSI'] = float(row.initial_surface_storage)
        sitedata['WAV'] = float(row.inital_water_availability)
        sitedata['SMLIM'] = float(row.smlim)
    else:
        raise RuntimeError("No rows found")

    return sitedata


class GridWeatherDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the GRID_WEATHER table in a CGMS database.

    :param metadata: SqlAlchemy metadata object providing DB access
    :param grid_no:  CGMS Grid ID
    :param startdate: Retrieve meteo data starting with startdate
        (datetime.date object)
    :param enddate: Retrieve meteo data up to and including enddate
        (datetime.date object)
    :param recalc_ET: Set to True to force calculation of reference
        ET values. Mostly useful when values have not been calculated
        in the CGMS database.
    :param use_cache: Set to False to ignore read/writing a cache file.

    Note that all meteodata is first retrieved from the DB and stored
    internally. Therefore, no DB connections are stored within the class
    instance. This makes that class instances can be pickled.

    """
    # default values for the Angstrom parameters in the sunshine duration model
    angstA = 0.29
    angstB = 0.49

    def __init__(self, DBconn, grid_no, start_date=None, end_date=None,
                 recalc_ET=False, use_cache=True):

        WeatherDataProvider.__init__(self)
        self.grid_no = int(grid_no)
        self.recalc_ET = recalc_ET
        self.use_cache = use_cache

        if not self._self_load_cache(self.grid_no) or self.use_cache is False:
            if start_date is None:
                start_date = dt.date(dt.MINYEAR, 1, 1)
            if end_date is None:
                end_date = dt.date(dt.MAXYEAR, 1, 1)
            self.start_date = self.check_keydate(start_date)
            self.end_date = self.check_keydate(end_date)
            self.timeinterval = (end_date - start_date).days + 1

            cursor = DBconn.cursor()

            # Get location info (lat/lon/elevation)
            self._fetch_location_from_db(cursor)

            # Retrieved meteo data
            self._fetch_grid_weather_from_db(cursor)

            # Description
            self.description = "Weather data derived for grid_no: %i" % self.grid_no

            # Save cache file
            if self.use_cache:
                fname = self._get_cache_filename(self.grid_no)
                self._dump(fname)

    def _get_cache_filename(self, grid_no):
        fname = "%s_grid_%i.cache" % (self.__class__.__name__, grid_no)
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, fname)
        return cache_filename

    def _self_load_cache(self, grid_no):
        """Checks if a cache file exists and tries to load it."""
        cache_fname = self._get_cache_filename(grid_no)
        if os.path.exists(cache_fname):
            r = os.stat(cache_fname)
            cache_file_date = dt.date.fromtimestamp(r.st_mtime)
            age = (dt.date.today() - cache_file_date).days
            if age < 1:
                try:
                    self._load(cache_fname)
                    return True
                except PCSEError:
                    pass
        return False

    def _fetch_location_from_db(self, cursor):
        """Retrieves latitude, longitude, elevation from 'grid' table and
        assigns them to self.latitude, self.longitude, self.elevation."""

        # Pull Latitude value for grid nr from database
        sql = "select latitude, longitude, altitude from grid where grid_no=?"
        r = cursor.execute(sql, (self.grid_no,))
        row = r.fetchone()
        if not row:
            raise PCSEError(f"Cannot find lat/lon for grid {self.grid_no}")

        self.latitude = float(row.latitude)
        self.longitude = float(row.longitude)
        self.elevation = float(row.altitude)

    def _fetch_grid_weather_from_db(self, cursor):
        """Retrieves the meteo data from table 'grid_weather'.
        """

        try:
            sql = "select * from grid_weather where grid_no=? and day>=? and day<=?"
            r = cursor.execute(sql, (self.grid_no, self.start_date, self.end_date))
            rows = r.fetchall()

            c = len(rows)
            if c < self.timeinterval:
                msg = "Only %i records selected from table 'grid_weather' "+\
                       "for grid %i, period %s -- %s."
                self.logger.warn(msg % (c, self.grid_no, self.start_date,
                                        self.end_date))

            meteopackager = self._make_WeatherDataContainer
            for row in rows:
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = meteopackager(row, t)
                self._store_WeatherDataContainer(wdc, DAY)
        except Exception as e:
            errstr = "Failure reading meteodata for day %s: %s" % (row.day, str(e))
            raise PCSEError(errstr)

    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""

        t.update({"TMAX": float(row.maximum_temperature),
                  "TMIN": float(row.minimum_temperature),
                  "VAP":  float(row.vapour_pressure),
                  "WIND": wind10to2(float(row.windspeed)),
                  "RAIN": float(row.rainfall)/10.,
                  "IRRAD": float(row.calculated_radiation)*1000.,
                  "SNOWDEPTH": safe_float(row.snow_depth)})

        if not self.recalc_ET:
            t.update({"E0":  float(row.e0)/10.,
                      "ES0": float(row.es0)/10.,
                      "ET0": float(row.et0)/10.})
        else:
            e0, es0, et0 = reference_ET(ANGSTA=self.angstA,
                                        ANGSTB=self.angstB, **t)
            t.update({"E0":  e0/10.,
                      "ES0": es0/10.,
                      "ET0": et0/10.})

        wdc = WeatherDataContainer(**t)

        return wdc
