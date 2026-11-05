# -*- coding: utf-8 -*-
# Copyright (c) 2026 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), November 2026
"""A weather data provider for reading AgERA5v2 data from the HTTP API provided by Wageningen-UR
at https://agera5.containers.wur.nl.
"""
import json

import requests
import pandas as pd

from pcse.base import WeatherDataContainer, WeatherDataProvider
from pcse.util import reference_ET, wind10to2, check_date
from pcse.exceptions import PCSEError


class AgERA5WeatherDataProvider(WeatherDataProvider):
    """WeatherDataProvider for using AgERA5 v2.0 data with PCSE.

    :param latitude: latitude [-180..180]to request weather data for
    :param longitude: longitude [-70..70] to request weather data for
    :param startdate: the start date for which to retrieve data [date, datetime or ISO date string]
    :param enddate: the end date for which to retrieve data [date, datetime or ISO date string]

    AgERA5 is a agronomically tailored version derived from the ECWMF ERA5 reanalysis.
    It provides 27 daily variables, although only 7 are used here as they are relevant for
    modelling with WOFOST or other PCSE models. These include:
    - temperature_max [C]
    - temperature_min [C]
    - temperature_avg [C]
    - vapourpressure [hPa]
    - windspeed [m2/sec 2m ]
    - precipitation [cm]
    - radiation [J/m2]
    - snowdepth [cm]

    AgERA5 is available from Climate Data Store here: https://doi.org/10.24381/f9c1f09d but this
    weather data provider uses the API available at https://agera5.containers.wur.nl

    This dataprovider does not provide caching for weather data.
    """
    no_conv = lambda x: x
    mm_to_cm = lambda x: x / 10.
    kJ_to_J = lambda x: x * 1000.
    to_date = lambda x: x.date()

    variable_renaming = [("temperature_max", "TMAX", no_conv),
                         ("temperature_min", "TMIN", no_conv),
                         ("temperature_avg",  "TEMP", no_conv),
                         ("vapourpressure", "VAP", no_conv),
                         ("windspeed", "WIND", wind10to2),
                         ("precipitation", "RAIN",mm_to_cm),
                         ("radiation", "IRRAD", no_conv),
                         ("snowdepth", "SNOWDEPTH", no_conv),
                         ("day", "DAY", to_date)]
    angstA = 0.25
    angstB = 0.45
    ETmodel = "PM"
    url = "https://agera5.containers.wur.nl/api/v1/get_agera5"


    def __init__(self, *, latitude=None, longitude=None, startdate=None, enddate=None):
        WeatherDataProvider.__init__(self)
        inputs = dict(latitude=latitude,
                      longitude=longitude,
                      startdate=check_date(startdate),
                      enddate=check_date(enddate)
                      )
        r = requests.get(self.url, params=inputs)
        data = json.loads(r.text)
        self.elevation = data["location_info"]["grid_agera5_elevation"]
        self.longitude = data["location_info"]["grid_agera5_longitude"]
        self.latitude = data["location_info"]["grid_agera5_latitude"]
        df = pd.DataFrame(data["weather_variables"])
        df["day"] = pd.to_datetime(df.day)
        df_new = self._rename_columns(df)
        for daily_weather in df_new.itertuples(index=False):
            self._make_WeatherDataContainer(daily_weather)

    def _rename_columns(self, df):
        """Renames columns and applies conversions for use in PCSE models
        """
        df_new = pd.DataFrame()
        for old_name, new_name, conversion in self.variable_renaming:
            df_new[new_name] = df[old_name].apply(conversion)

        return df_new

    def _make_WeatherDataContainer(self, daily_weather):
        """Converts a record of daily weather data into a WeatherDataContainer and stores it.
        """
        t = daily_weather._asdict()
        t.update({"LAT": self.latitude, "LON": self.longitude, "ELEV": self.elevation})

        # Reference evapotranspiration in mm/day
        try:
            E0, ES0, ET0 = reference_ET(ANGSTA=self.angstA, ANGSTB=self.angstB, ETMODEL=self.ETmodel, **t)
        except ValueError as e:
            msg = (f"Failed to calculate reference ET values on {daily_weather.DAY}. " +
                   f"With input values:\n {t}\n Due to error: {e}")
            raise PCSEError(msg)

        # update record with ET values value convert to cm/day
        t.update({"E0": E0 / 10., "ES0": ES0 / 10., "ET0": ET0 / 10.})

        # Build and store weather data container from dict 't' for this day
        wdc = WeatherDataContainer(**t)
        self._store_WeatherDataContainer(wdc, daily_weather.DAY)


class AgERA5CHIRPSWeatherDataProvider(AgERA5WeatherDataProvider):
    """WeatherDataProvider for using AgERA5 v2.0 data. Rainfall from AgERA5 is replaced with
    rainfall from CHIRPS v3. Particularly in the tropics CHIRPS has better rainfall estimates compared to AgERA5.

    :param latitude: latitude [-180..180]to request weather data for
    :param longitude: longitude [-70..70] to request weather data for
    :param startdate: the start date for which to retrieve data [date, datetime or ISO date string]
    :param enddate: the end date for which to retrieve data [date, datetime or ISO date string]

    AgERA5 is available from Climate Data Store here: https://doi.org/10.24381/f9c1f09d
    CHIRPS v3 is available from UCSB: https://www.chc.ucsb.edu/data/chirps3

    This weather data provider uses the API available at https://agera5.containers.wur.nl

    This dataprovider does not provide caching for weather data.
    """
    url = "https://agera5.containers.wur.nl/api/v1/get_agera5_chirps"
