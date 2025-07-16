# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
"""
Data providers for weather, agromanagement, soil, crop and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 8 database schema.
"""
import os
import datetime as dt

from sqlalchemy import MetaData, select, Table, and_
import yaml

from ...util import check_date, wind10to2, reference_ET, safe_float
from ... import settings
from ... import exceptions as exc
from ..cgms12.data_providers import fetch_crop_name
from ...base import WeatherDataProvider, WeatherDataContainer

#----------------------------------------------------------------------------
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

    def __init__(self, engine, grid_no, start_date=None, end_date=None,
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

            meta = MetaData()
            meta.reflect(bind=engine)

            # Get location info (lat/lon/elevation)
            self._fetch_location_from_db(engine, meta)

            # Retrieved meteo data
            self._fetch_grid_weather_from_db(engine, meta)

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
                except exc.PCSEError:
                    pass
        return False

    #---------------------------------------------------------------------------
    def _fetch_location_from_db(self, engine, meta):
        """Retrieves latitude, longitude, elevation from 'grid' table and
        assigns them to self.latitude, self.longitude, self.elevation."""

        # Pull Latitude value for grid nr from database

        table_grid = meta.tables['grid']
        s = select(table_grid.c.latitude, table_grid.c.longitude, table_grid.c.altitude)
        s = s.where(table_grid.c.grid_no==self.grid_no)
        with engine.connect() as con:
            cursor = con.execute(s)
            row = cursor.fetchone()

        if row is None:
            msg = "Failed deriving location info for grid %s: %s" % (self.grid_no, e)
            raise exc.PCSEError(msg)

        self.latitude = row.latitude
        self.longitude = row.longitude
        self.elevation = row.altitude

        msg = "Succesfully retrieved location information from 'grid' table "+\
              "for grid %s"
        self.logger.info(msg % self.grid_no)

    def _fetch_grid_weather_from_db(self, engine, metadata):
        """Retrieves the meteo data from table 'grid_weather'.
        """

        try:
            table_gw = Table('grid_weather', metadata, autoload=True)
            s = select(table_gw)
            s = s.where(and_(table_gw.c.grid_no==self.grid_no,
                             table_gw.c.day>=self.start_date,
                             table_gw.c.day<=self.end_date)
                       )
            with engine.connect() as DBconn:
                cursor = DBconn.execute(s)
                rows = cursor.fetchall()

            c = len(rows)
            if c < self.timeinterval:
                msg =  "Only %i records selected from table 'grid_weather' "+\
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
            raise exc.PCSEError(errstr)

        msg = ("Successfully retrieved weather data from 'grid_weather' table "
               "for grid %s between %s and %s")
        self.logger.info(msg % (self.grid_no, self.start_date, self.end_date))

    #---------------------------------------------------------------------------
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



class AgroManagementDataProvider(list):
    """Class for providing agromanagement data from the CROP_CALENDAR table in a CGMS8 database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the grid_no column in the table
    :param crop_no: Integer id of crop, maps to the crop_no column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
           The campaign year refers to the year of the crop start. Thus for crops
           crossing calendar years, the start_date can be in the previous year as the
           harvest.

    Note that by default the campaign_start_date is set equal to the crop_start_date which
    means that the simulation starts when the crop starts. In some cases this is undesirable
    and an earlier start date should be used. In that case use the `set_campaign_start_date(date)`
    to update the campaign_start_date.
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

    def __init__(self, engine, grid_no, crop_no, campaign_year):
        list.__init__(self)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.crop_name = fetch_crop_name(engine, self.crop_no)
        self.db_resource = str(engine)[7:-1]

        metadata = MetaData(engine)
        table_cc = Table("crop_calendar", metadata, autoload=True)

        r = select([table_cc], and_(table_cc.c.grid_no == self.grid_no,
                                    table_cc.c.crop_no == self.crop_no,
                                    table_cc.c.year == self.campaign_year)).execute()
        row = r.fetchone()
        r.close()
        if row is None:
            msg = "Failed deriving crop calendar for grid_no %s, crop_no %s " % (grid_no, crop_no)
            raise exc.PCSEError(msg)

        # Determine the start date.
        year = int(row.year)
        month = int(row.start_month1)
        day = int(row.start_monthday1)
        self.crop_start_date = check_date(dt.date(year, month, day))
        self.campaign_start_date = self.crop_start_date

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        cgms_start_type = str(row.start_type).strip()
        if cgms_start_type == "fixed_sowing":
            self.crop_start_type = "sowing"
        elif cgms_start_type == "fixed_emergence":
            self.crop_start_type = "emergence"
        else:
            msg = "Unsupported START_TYPE in CROP_CALENDAR table: %s" % row.start_type
            raise exc.PCSEError(msg)

        # Determine maximum duration of the crop
        self.max_duration = int(row.max_duration)

        # Determine crop end date/type and the end of the campaign
        self.crop_end_type = str(row.end_type).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % row.end_type)
            raise exc.PCSEError(msg)

        # Note that we end one day to the campaign end date in order to avoid that the
        # crop_end_date and campaign_end_date fall on the same date.
        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
            self.campaign_end_date = self.crop_start_date + dt.timedelta(days=self.max_duration+1)
        else:
            month = int(row.end_month)
            day = int(row.end_monthday)
            self.crop_end_date = dt.date(year, month, day)
            if self.crop_end_date <= self.crop_start_date:
                self.crop_end_date = dt.date(year+1, month, day)
            self.campaign_end_date = self.crop_end_date + dt.timedelta(days=1)

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
                                                     duration=self.max_duration,
                                                     campaign_end_date=self.campaign_end_date
                                                     )
        return input

    def _parse_yaml(self, input):
        """Parses the input YAML string and assigns to self"""
        try:
            items = yaml.safe_load(input)
        except yaml.YAMLError as e:
            msg = "Failed parsing agromanagement string %s: %s" % (input, e)
            raise exc.PCSEError(msg)
        del self[:]
        self.extend(items)

    def set_campaign_start_date(self, start_date):
        """Updates the value for the campaign_start_date.

        This is useful only when the INITIAL_SOIL_WATER table in CGMS8 defines a different
        campaign start date which should be used instead.
        """
        self.campaign_start_date = check_date(start_date)
        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)
