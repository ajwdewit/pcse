# -*- coding: utf-8 -*-
# CC-BY-SA, Alterra, Wageningen-UR
# Jappe Franke (jappe.franke@wur.nl), Allard de Wit (allard.dewit@wur.nl), June 2016
"""A weather data provider reading its data from GeoBIS mysql database.
Data providers are compatible with a CGMS 8 database schema.
"""
import datetime as dt

import pymysql
from pandas.core.frame import DataFrame

from .. import config
from .. import exceptions as exc
from ..base_classes import WeatherDataProvider, WeatherDataContainer
from ..util import wind10to2



#----------------------------------------------------------------------------
class GridWeatherDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the GRID_WEATHER table in a CGMS(like) database.

    Note that all meteodata is first retrieved from the DB and stored
    internally. Therefore, no DB connections are stored within the class
    instance. This makes that class instances can be pickled.
    """

    def __init__(self, lat, lon):

        WeatherDataProvider.__init__(self)

        start_date = dt.date(dt.MINYEAR, 1, 1)
        end_date = dt.date(dt.MAXYEAR, 1, 1)
        self.start_date = self.check_keydate(start_date)
        self.end_date = self.check_keydate(end_date)
        self.timeinterval = (end_date - start_date).days + 1

        self.latitude = lat
        self.longitude = lon
        self.elevation = None
        
        #set DB connect, we use pymysql here for calling stored procedures
        self.connection = pymysql.connect(host= config.myhost, user= config.myuser, password=config.mypwd,
                                          db= config.mydb, cursorclass=pymysql.cursors.DictCursor)

        # Get location info (lat/lon/elevation)
        self._fetch_location_from_db()

        # Retrieved meteo data
        self._fetch_grid_weather_from_db()

        self.connection.close()

        self.description = "Weather data provider for ISIDORE."
    #---------------------------------------------------------------------------
    
    def _fetch_location_from_db(self):
        """Retrieves grid_no from 'grid' table via  a stored procedure and
        assigns it to self.grid_no"""
        # Pull grid nr from database

        try:
            cur = self.connection.cursor()
            cur.execute("call get_grid(%s,%s)" % (self.latitude,self.longitude));
            row = cur.fetchall()
            cur.close()
            if not row:
                raise Exception()
            self.grid_no = row[0]['grid_no']
        except Exception as e:
            msg = "Failed deriving grid info for lat %s, lon %s" % (self.latitude,self.longitude)
            raise exc.PCSEError(msg)

        try:
            cur = self.connection.cursor()
            cur.execute("select altitude from grid where grid_no=%s" % self.grid_no)
            row = cur.fetchone()
            self.elevation = row['altitude']
        except Exception as e:
            msg = "Failed deriving altitude info for grid %s" % (self.grid_no)
            raise exc.PCSEError(msg)

        msg = "Succesfully retrieved location information for 'grid' table "+\
              "for grid %s"
        #self.logger.info(msg % self.grid_no)
    #---------------------------------------------------------------------------
    def _fetch_grid_weather_from_db(self):
        """Retrieves the meteo data from stored procedure 'grid_weather'.
        """

        try:
            cur = self.connection.cursor()
            cur.execute("call get_grid_weather(%s)"%(self.grid_no));
            rows = DataFrame(cur.fetchall())
            cur.close()
            c = len(rows)
            if c < self.timeinterval:
                msg =  "Only %i records selected from table 'grid_weather' "+\
                       "for grid %i, period %s -- %s."
                #self.logger.warn(msg % (c, self.grid_no, self.start_date, self.end_date))
    
            meteopackager = self._make_WeatherDataContainer
            for row in rows.itertuples():
                if row.day is None:
                    continue
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = meteopackager(row, t)
                self._store_WeatherDataContainer(wdc, DAY)
        except Exception as e:
            errstr = "Failure reading meteodata: " + str(e)
            raise exc.PCSEError(errstr)

        msg = ("Successfully retrieved weather data from 'grid_weather' table "
               "for grid %s between %s and %s")
        #self.logger.info(msg % (self.grid_no, self.start_date, self.end_date))
        
    #---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""

        t.update({"TMAX": float(row.maximum_temperature),
                  "TMIN": float(row.minimum_temperature),
                  "VAP":  float(row.vapour_pressure),
                  "WIND": wind10to2(float(row.windspeed)),
                  "RAIN": float(row.rainfall)/10.,
                  "E0":  float(row.e0)/10.,
                  "ES0": float(row.es0)/10.,
                  "ET0": float(row.et0)/10.,
                  "IRRAD": float(row.calculated_radiation)})
        wdc = WeatherDataContainer(**t)

        return wdc
    #---------------------------------------------------------------------------


    def export(self):
        """Exports the contents of the WeatherDataProvider as a list of dictionaries.

        The results from export can be directly converted to a Pandas dataframe
        which is convenient for plotting or analyses.
        """
        if self.supports_ensembles:
            # We have to include the member_id in each dict with weather data
            pass
        else:
            weather_data = []
            days = sorted([r[0] for r in self.store.keys()])
            for day in days:
                wdc = self(day)
                r = {key: getattr(wdc, key) for key in wdc.__slots__ if hasattr(wdc, key)}
                weather_data.append(r)
        return weather_data


if __name__ == "__main__":
    latitude = 25.733
    longitude = 91.2333
    wdp= GridWeatherDataProvider(latitude, longitude)
    print "Done"
    pass
