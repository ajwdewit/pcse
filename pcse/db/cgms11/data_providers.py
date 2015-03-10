"""Data providers for weather, soil, timer and site compatible with a
CGMS 11 compatible database.
"""
import sys, os
import datetime

from sqlalchemy import create_engine, MetaData, select, Table, and_, join

from ...util import wind10to2, safe_float
from ... import exceptions as exc
from ...base_classes import WeatherDataContainer, WeatherDataProvider


class WeatherObsGridDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the WEATHER_OBS_GRID table in a CGMS11
    compatible database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no:  Grid ID of PyWofost run
    :param start_date: Retrieve meteo data starting with start_date
        (datetime.date object)
    :param end_date: Retrieve meteo data up to and including end_date
        (datetime.date object)

    Note that all meteodata is first retrieved from the DB and stored
    internally. Therefore, no DB connections are stored within the class
    instance. This makes that class instances can be pickled.

    If start_date and end_date are not provided then the entire time-series
    for the grid is retrieved.
    """

    def __init__(self, engine, grid_no, start_date=None, end_date=None):

        WeatherDataProvider.__init__(self)

        self.grid_no = grid_no
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

        metadata = MetaData(engine)
        # Get location info (lat/lon/elevation)
        self._fetch_location_from_db(metadata)
        # Retrieved meteo data
        self._fetch_weather_from_db(metadata)

        # Provide a description that is shown when doing a print()
        line1 = "Weather data retrieved from CGMS 11 db %s" % str(engine)[7:-1]
        line2 = "for grid_no: %s" % self.grid_no
        self.description = [line1, line2]

    #---------------------------------------------------------------------------
    def _fetch_location_from_db(self, metadata):
        """Retrieves latitude, longitude, elevation from 'grid' table and
        assigns them to self.latitude, self.longitude, self.elevation."""

        # Pull Latitude value for grid nr from database

        table_grid = Table('grid', metadata, autoload=True)
        r = select([table_grid.c.latitude, table_grid.c.longitude,
                    table_grid.c.altitude],
                   table_grid.c.grid_no == self.grid_no).execute()
        row = r.fetchone()
        r.close()
        if row is None:
            msg = "Failed deriving location info for grid %s" % self.grid_no
            raise exc.PCSEError(msg)

        self.latitude = float(row.latitude)
        self.longitude = float(row.longitude)
        self.elevation = float(row.altitude)

        msg = ("Succesfully retrieved location information from 'grid' table "
               "for grid %s")
        self.logger.info(msg, self.grid_no)

    #---------------------------------------------------------------------------
    def _fetch_weather_from_db(self, metadata):
        """Retrieves the meteo data from table 'grid_weather'.
        """

        try:
            # if start_date/end_date are None, define a date in the far past/future
            start_date = self.start_date if self.start_date is not None else datetime.date(1, 1, 1)
            end_date = self.end_date if self.end_date is not None else datetime.date(9999, 1, 1)
            table_db = Table('weather_obs_grid', metadata, autoload=True)
            r = select([table_db], and_(table_db.c.grid_no == self.grid_no,
                                        table_db.c.day >= start_date,
                                        table_db.c.day <= end_date)
                       ).execute()
            rows = r.fetchall()

            c = len(rows)
            if self.time_interval is not None:
                if c < self.time_interval:
                    msg = ("Only %i records selected from table 'grid_weather' "
                           "for grid %i, period %s -- %s.")
                    self.logger.warn(msg, c, self.grid_no, self.start_date,
                                     self.end_date)

            for row in rows:
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = self._make_WeatherDataContainer(row, t)
                self._store_WeatherDataContainer(wdc, DAY)
        except Exception as e:
            msg = "Failure reading meteodata for grid %s "
            self.logger.exception(msg, self.grid_no)
            raise exc.PCSEError(msg, self.grid_no)

        msg = ("Successfully retrieved weather data from 'grid_weather' table "
               "for grid %s between %s and %s")
        self.logger.info(msg, self.grid_no, self.start_date, self.end_date)

    #---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""

        t.update({"TMAX": float(row.temperature_max),
                  "TMIN": float(row.temperature_min),
                  "TEMP": float(row.temperature_avg),
                  "VAP":  float(row.vapourpressure),
                  "WIND": wind10to2(float(row.windspeed)),
                  "RAIN": float(row.precipitation)/10.,
                  "E0":  float(row.e0)/10.,
                  "ES0": float(row.es0)/10.,
                  "ET0": float(row.et0)/10.,
                  "IRRAD": float(row.radiation)*1000.,
                  "SNOWDEPTH": safe_float(row.snowdepth)})
        wdc = WeatherDataContainer(**t)
        return wdc
