# -*- coding: utf-8 -*-
# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2016

"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 14 database schema.
"""

import datetime

from sqlalchemy import MetaData, select, Table, and_
# from tabulate import tabulate
import numpy as np

from ...util import wind10to2, safe_float, check_date, reference_ET
from ... import exceptions as exc
from ...base_classes import WeatherDataContainer, WeatherDataProvider

def fetch_crop_name(engine, idcrop):
    """Retrieves the name of the crop from the CROP table for given idcrop.

    :param engine: SqlAlchemy engine object providing DB access
    :param idcrop: Integer crop ID, maps to the idcrop column in the table
    """
    result = "Unknown"
    metadata = MetaData(engine)
    table_crop = Table("global_crops", metadata, autoload=True)
    sc = select([table_crop], table_crop.columns.idcrop == idcrop).execute()
    row = sc.fetchone()
    sc.close()
    if row is None:
        msg = "Failed deriving crop name from view GLOBAL_CROPS for idcrop %s" % idcrop
        raise exc.PCSEError(msg)
    result = row.descriptor
    return result

class STU_Suitability(set):
    """Returns a set() of suitable STU's for given idcrop.

    :param engine: SqlAlchemy engine object providing DB access
    :param idcrop: Integer crop ID, maps to the idcrop column in the table
    """
    def __init__(self, engine, idcrop):
        # Initialise
        self.idcrop = int(idcrop)
        self.crop_name = fetch_crop_name(engine, idcrop)
        metadata = MetaData(engine)

        # Retrieve suitable STU's        
        table_link = Table("link_crop_stu", metadata, autoload=True)
        sc = select([table_link], table_link.columns.idcrop == idcrop).execute()
        rows = sc.fetchall()
        sc.close()
        if rows is None or len(rows) == 0:
            msg = "No suitable soil type unit found for idcrop=%s" % idcrop
            raise exc.PCSEError(msg)
        set.__init__(self, [int(row.idstu) for row in rows])

class WeatherObsGridDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the WEATHER_ERA_GRID table in a CGMS14
    compatible database - or from another table when specified

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid:  Grid number (int) to retrieve data for
    :param start_date: Retrieve meteo data starting with start_date
        (datetime.date object)
    :param end_date: Retrieve meteo data up to and including end_date
        (datetime.date object)
    :param recalc_ET: Set to True to force calculation of reference
        ET values. Mostly useful when values have not been calculated
         in the CGMS database.
    :param table_name: 
    """
        # default values for the Angstrom parameters in the sunshine duration model
    angstA = 0.18
    angstB = 0.55
   
    def __init__(self, engine, idgrid, start_date=None, end_date=None, 
            recalc_ET=False, table_name='weather_era_grid'):        # Initialise
        WeatherDataProvider.__init__(self)
        self.idgrid = idgrid
        self.recalc_ET = recalc_ET
        self.table_name = table_name
        metadata = MetaData(engine)
        
        # Check the start and end dates and assign when ok
        try:
            self.start_date = self.check_keydate(start_date)
        except KeyError:
            self.start_date = None
        try:
            self.end_date = self.check_keydate(end_date)
        except KeyError:
            self.end_date = None
        try:
            self.time_interval = (end_date - start_date).days + 1
        except TypeError:
            self.time_interval = None
        
        # Get location info (lat/lon/elevation)
        self._fetch_location_from_db(metadata)
        
        # Retrieve the meteo data
        self._fetch_weather_from_db(metadata)
        
        # Provide a description that is shown when doing a print()
        line1 = "Weather data retrieved from CGMS 14 db %s" % str(engine)[7:-1]
        line2 = "for idgrid: %s" % self.idgrid
        self.description = [line1, line2]
        
    #---------------------------------------------------------------------------
    def _fetch_location_from_db(self, metadata):
        """Retrieves latitude, longitude, elevation from "grid" table and
        assigns them to self.latitude, self.longitude, self.elevation."""
        tg = Table("grids", metadata, autoload=True)
        sc = select([tg.c.latitude, tg.c.longitude, tg.c.altitude],
                tg.c.idgrid == self.idgrid).execute()
        row = sc.fetchone()
        sc.close()
        if row is None:
            msg = "Failed deriving location info for grid %s" % self.idgrid
            raise exc.PCSEError(msg)
        
        # Use the resulta
        self.latitude = float(row.latitude)
        self.longitude = float(row.longitude)
        self.elevation = float(row.altitude)

        # Report success
        msg = ("Succesfully retrieved location information from 'grids' table for grid %s")
        self.logger.info(msg, self.idgrid)
        
    #---------------------------------------------------------------------------
    def _fetch_weather_from_db(self, metadata):
        """Retrieves the meteo data from table "grid_weather".
        """
        try:
            # Check input - if start_date/end_date are None, define a date in the far past/future
            start_date = self.start_date if self.start_date is not None else datetime.date(1, 1, 1)
            end_date = self.end_date if self.end_date is not None else datetime.date(9999, 1, 1)
            
            # Get hold of weather table, statement, sqlalchemy cursor and finally of the rows
            tw = Table(self.table_name, metadata, autoload=True)
            sm = select([tw], 
                and_(tw.c.idgrid == self.idgrid, tw.c.day >= start_date, tw.c.day <= end_date))
            sc = sm.execute()
            rows = sc.fetchall()
        
            count = len(rows)
            if self.time_interval is not None:
                if count < self.time_interval:
                    msg = ("Only %i records selected from table 'WEATHER_OBS_GRID' "
                           "for grid %i, period %s -- %s.")
                    self.logger.warn(msg, count, self.idgrid, self.start_date,
                                     self.end_date)

            for row in rows:
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = self._make_WeatherDataContainer(row, t)
                self._store_WeatherDataContainer(wdc, DAY)
        
        except Exception:
            msg = "Failure reading meteodata for grid %s "
            self.logger.exception(msg, self.idgrid)
            raise exc.PCSEError(msg, self.idgrid)

        # Report success
        msg = ("Successfully retrieved weather data from " + self.table_name + " table "
               "for grid %s between %s and %s")
        self.logger.info(msg, self.idgrid, self.start_date, self.end_date)
    
    #---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""
        result = None

        t.update({"TMAX": float(row.temperature_max),
                  "TMIN": float(row.temperature_min),
                  "TEMP": float(row.temperature_avg),
                  "VAP":  float(row.vapourpressure),
                  "WIND": wind10to2(float(row.windspeed)),
                  "RAIN": float(row.precipitation)/10.,
                  "IRRAD": float(row.radiation)*1000.,
                  "SNOWDEPTH": safe_float(row.snowdepth)})

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

        result = WeatherDataContainer(**t)
        return result
    
class TimerDataProvider(dict):
    """Class for providing timerdata from the CROP_CALENDAR table in a CGMS14 database.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid: Integer grid ID, maps to the idgrid column in the table
    :param idcrop: Integer crop ID, maps to the idcrop column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """
    def __init__(self, engine, idgrid, idcrop, campaign_year): 
        # Initialise  
        self.idgrid = int(idgrid)
        self.idcrop = int(idcrop) 
        self.campaign_year = int(campaign_year)
        metadata = MetaData(engine)
    
    def set_START_DATE(self, start_date):
        """Updates the value for START_DATE in TimerDataProvider
        """
        pass
    
    def __str__(self):
        result = ""
        
        return result

class SoilDataProviderSingleLayer(dict):
    """Class for providing soil data from the ROOING_DEPTH AND
    SOIL_PHYSICAL_GROUP tableS in a CGMS14 database. This
    applies to the single layered soil only.

    :param engine: SqlAlchemy engine object providing DB access
    :param idstu: Integer stu no, maps to the idstu column in the table
                   SOIL_TYPOLOGIC_UNIT for providing soil data

    Note that the value of the parameter SMLIM (Initial maximum
    moisture content in initial rooting depth zone) is set
    to field capacity (SMFCF)
    """
    soil_parameters = [("CRAIRC", "CRITICAL_AIR_CONTENT"),
        ("K0", "HYDR_CONDUCT_SATUR"),
        ("SOPE", "MAX_PERCOL_ROOT_ZONE"),
        ("KSUB", "MAX_PERCOL_SUBSOIL"),
        ("SMFCF", "SOIL_MOISTURE_CONTENT_FC"),
        ("SM0", "SOIL_MOISTURE_CONTENT_SAT"),
        ("SMW", "SOIL_MOISTURE_CONTENT_WP")]

    def __init__(self, engine, idstu):
        dict.__init__(self)
        metadata = MetaData(engine)
    
    def _get_from_STU(self, metadata, idstu):
        """Retrieves the soil parameters for the given soil typologic unit
        (idstu) from the tables SOIL_PHYSICAL_GROUP and ROOTING_DEPTH.
        """
        rd_class = 0
        soil_group_no = 0
        table_stu = Table("soil_stu", metadata, autoload=True)
        
        return (rd_class, soil_group_no)
    
    def _get_rooting_depth(self, metadata, rd_class):
        """Gets the rooting depth from the table ROOTING_DEPTH and
         stores into self[] directly under parameter name 'RDMSOL'.

        :param metadata: An SQLAlchemy Metadata object
        :param rd_class: The rooting depth class (integer)
        """
        pass
    
    def _get_soil_hydraulic_parameters(self, metadata, spg_no):
        """Retrieves the soil hydraulic parameters and stores into self[] directly.

        :param metadata: An SQLAlchemy Metadata object
        :param spg_no: the soil physical group number (integer)
        :return: None
        """
        pass
    
class SoilDataIterator(list):
    """Class for iterating over the different soils in a CGMS grid.

    Instances of this class behave like a list, allowing to iterate
    over the soils in a CGMS grid. An example::

    >>> soil_iterator = SoilDataIterator(engine, idgrid=15060)
    >>> print(soildata)
    Soil data for idgrid=15060 derived from oracle+cx_oracle://cgms12eu:***@eurdas.world
      smu_no=9050131, area=625000000, idstu=9000282 covering 50% of smu.
        Soil parameters {'SMLIM': 0.312, 'SMFCF': 0.312, 'SMW': 0.152, 'CRAIRC': 0.06,
                         'KSUB': 10.0, 'RDMSOL': 10.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.439}
      smu_no=9050131, area=625000000, idstu=9000283 covering 50% of smu.
        Soil parameters {'SMLIM': 0.28325, 'SMFCF': 0.28325, 'SMW': 0.12325, 'CRAIRC': 0.06,
                         'KSUB': 10.0, 'RDMSOL': 40.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.42075}
    >>> for smu_no, area, idstu, percentage, soil_par in soildata:
    ...     print(smu_no, area, idstu, percentage)
    ...
    (9050131, 625000000, 9000282, 50)
    (9050131, 625000000, 9000283, 50)
    """
    
    def __init__(self, engine, idgrid):
        list.__init__(self)
        self.idgrid = int(idgrid)

    def _get_SMU_from_EMU(self, metadata, idgrid):
        """Retrieves the relevant SMU for given idgrid from
        table EMU.
        """
        result = None
        
        if result is None:
            msg = ("No soil mapping units (SMU) found for idgrid=%i in table EMU" % idgrid)
            raise exc.PCSEError(msg)
        return result
    
    def _get_STU_from_SMU(self, metadata, smu_no):
        """Retrieves the relevant STU for given SMU_NO from table
        SOIL_ASSOCIATION_COMPOSITION
        """
        result = None
        
        if result is None:
            msg = "No soil typologic units (STU) found for smu_no=%i" % smu_no
            raise exc.PCSEError(msg)
        return result
    
    def __str__(self):
        result = ""
        
        return result
    
class CropDataProvider(dict):
    """Retrieves the crop parameters for the given idgrid, idcrop and year
    from the tables CROP_CALENDAR, CROP_PARAMETER_VALUE and VARIETY_PARAMETER_VALUE.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid: Integer grid ID, maps to the idgrid column in the table
    :param idcrop: Integer crop ID, maps to the idcrop column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """
    # Define single and tabular crop parameter values
    parameter_codes_single = ("CFET", "CVL", "CVO", "CVR", "CVS", "DEPNR", "DLC",
        "DLO", "DVSEND", "EFF", "IAIRDU", "IDSL", "KDIF",
        "LAIEM", "PERDL", "Q10", "RDI", "RDMCR", "RGRLAI",
        "RML", "RMO", "RMR", "RMS", "RRI", "SPA", "SPAN", "SSA",
        "TBASE", "TBASEM", "TDWI", "TEFFMX", "TSUM1", "TSUM2",
        "TSUMEM")
    parameter_codes_tabular = ("AMAXTB", "DTSMTB", "FLTB", "FOTB", "FRTB", "FSTB",
        "RDRRTB", "RDRSTB", "RFSETB", "SLATB", "TMNFTB",
         "TMPFTB")
    # Some parameters have to be converted from a single to a tabular form
    single2tabular = {"SSA": ("SSATB", [0., None, 2.0, None]),
        "KDIF": ("KDIFTB", [0., None, 2.0, None]),
        "EFF": ("EFFTB", [0., None, 40., None])}
    # Default values for additional parameters not defined in CGMS
    parameters_additional = {"DVSI": 0.0, "IOX": 0}

    def __init__(self, engine, idgrid, idcrop, campaign_year):
        dict.__init__(self)
        self.idgrid = int(idgrid)
        self.idcrop = int(idcrop)
        self.campaign_year = int(campaign_year)
        
    def _fetch_crop_parameter_values(self, metadata, idcrop):
        """Derived the crop parameter values from the CROP_PARAMETER_VALUE
        table for given idcrop and add directly to dict self[]..
        """
        pass
    
    def _fetch_variety_parameter_values(self, metadata, idcrop, variety_no):
        """Derived the crop parameter values from the VARIETY_PARAMETER_VALUE
        table for given idcrop & variety_on and add directly to dict self[].
        """
        pass
    
    def _convert_single2tabular(self, parameter_code, pvalue):
        """Converts the single parameter into a tabular parameter.
        """
        tabular_parameter_code = ""
        tabular_values = []
        
        return tabular_parameter_code, tabular_values
    
    def __str__(self):
        result = ""
        
        return result
	
class SiteDataProvider(dict):
    """Provides the site data from the tables INITIAL_SOIL_WATER and SITE.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid:  Grid number (int)
    :param idcrop: Crop number (int)
    :param campaign_year: Campaign year (int)
    :param idstu: soil typologic unit number (int)

    Note that the parameter SSI (Initial surface storage) is
    set to zero

    """
    
    def __init__(self, engine, idgrid, idcrop, campaign_year, idstu):
        # Initialise
        dict.__init__(self)
        self.idgrid = int(idgrid)
        self.idcrop = int(idcrop)
        self.campaign_year = int(campaign_year)
        self.idstu = int(idstu)
        
        # Start using the engine
        self.db_resource = str(engine)[7:-1]
        self.crop_name = fetch_crop_name(engine, self.idcrop)
        metadata = MetaData(engine)
        
        table_crop_agg = Table('crop_aggregations', metadata, autoload=True)
        sc = select([table_crop_agg], table_crop_agg.c.idcrop == self.idcrop).execute()
        row = sc.fetchone()
        sc.close()
        if row == None:
            msg = ("Failed retrieving site data for grid_no=%s, crop_no=%s, campaign_year=%s, "
                "stu_no=%s" % (self.idgrid, self.idcrop, self.campaign_year, self.idstu))
            raise exc.PCSEError(msg)
        idcrop_parametrization = row.idcrop_parametrization
        t = Table('soil_initial_water', metadata, autoload=True)
        sm = select([t], 
            and_(t.c.idgrid == self.idgrid, t.c.idcrop_parametrization == idcrop_parametrization,
                 t.c.year == self.campaign_year, t.c.idstu == self.idstu))
        sc = sm.execute()
        
        # TODO: what can we learn from table soil_initial_water???
        
        # TODO: where can we find more site information?
        
    def __str__(self):
        result = ("Site parameter values for grid_no=%s, crop_no=%s (%s), stu_no=%s, "
            "campaign_year=%i derived from %s\n" % (self.grid_no, self.crop_no,
            self.crop_name, self.stu_no, self.campaign_year, self.db_resource))
        result += "    %s" % dict.__str__(self)
        return result
    