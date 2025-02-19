# -*- coding: utf-8 -*-
# Artificial Intelligence Group, WUR
# Hilmy Baja (hilmy.baja@wur.nl), February 2025
# A lot of code borrowed from nasapower.py by Allard de Wit

import os
import datetime
import time

from typing import Union

import requests

import pandas as pd
import numpy as np

from pcse.base import WeatherDataProvider, WeatherDataContainer
from pcse.util import reference_ET, wind10to2, check_angstromAB
from pcse.exceptions import PCSEError
from pcse.settings import settings


class OpenMeteoWeatherDataProvider(WeatherDataProvider):
    """
    A weather provider that uses the Open Meteo weather API.
    This object only needs a location (latitude and longitude)
    at initialization.
    There are two important parameters when constructing the object:
    :openmeteo_model and :forecast.

    The class variables list possible models to use, either for forecasts
    or historical data.

    To utilize a specific model, call it with the appropriate key argument.
    Be aware that there might be some nuances with using certain models.
    This hasn't been tested thoroughly, so there might be some issues with the starting
    date. Please provide an argument for the :start_date parameter if you find any issues.
    More info for each model is documented here: https://open-meteo.com/en/docs

    If you don't specify a model, the Open Meteo API will automatically
    choose the best model for your chosen location.
    """

    # Some class Variables
    HTTP_OK = 200
    angstA = 0.29
    angstB = 0.49

    #  List of forecast and historical weather models for OpenMeteo
    #  Comments show coverage and spatial resolution
    dict_forecast_models = {
        "best_match": datetime.date(2023, 1, 1),  # global, 0.25deg
        "arpae_cosmo_5m": datetime.date(2024, 2, 2),  # europe, 5m
        "bom_access_global": datetime.date(2024, 1, 19),  # global, 0.15deg
        "gem_seamless": datetime.date(2022, 11, 24),  # global, 0.15deg
        "jma_gsm": datetime.date(2016, 1, 1),  # global, 0.5deg
        "icon_seamless": datetime.date(2022, 11, 25),  # global, 11km
        "ecmwf_ifs025": datetime.date(2024, 2, 4),  # global, 0.25deg
        "knmi_seamless": datetime.date(2024, 7, 2),  # europe, 2.5km
        "meteofrance_seamless": datetime.date(2024, 1, 3),  # global 0.25deg
        "gfs_seamless": datetime.date(2021, 3, 24),  # global 0.11deg
        "ukmo_seamless": datetime.date(2022, 3, 2),  # global, 0.09deg/10km
    }

    dict_historical_models = {
        "best_match": datetime.date(1941, 1, 1), # global, 0.25deg
        "era5": datetime.date(1941, 1, 1),  # global, 0.25deg
        "era5_land": datetime.date(1951, 1, 1),  # global, 0.1deg
        "ecmwf_ifs": datetime.date(2017, 1, 1),  # global, 9km
        "cerra": datetime.date(1986, 1, 1),  # global, 5km
    }

    delay_historical_models = {
        "best_match": 10,
        "era5": 6,  # global, 0.25deg
        "era5_land": 6,  # global, 0.1deg
        "ecmwf_ifs": 3,  # global, 9km
        "cerra": 1,
    }

    def __init__(
        self,
        latitude: float,
        longitude: float,
        timezone: str = "UTC",
        openmeteo_model: str = "best_match",
        start_date: Union[str, datetime.date] = None,
        ETmodel: str = "PM",
        forecast: bool = False,
        force_update: bool = False,
    ):
        WeatherDataProvider.__init__(self)

        self.model = openmeteo_model
        self.ETmodel = ETmodel
        self.start_date = start_date
        self.is_forecast = forecast

        # update start date if using a specific model
        self._check_start_date()

        if latitude < -90 or latitude > 90:
            msg = "Latitude should be between -90 and 90 degrees."
            raise ValueError(msg)
        if longitude < -180 or longitude > 180:
            msg = "Longitude should be between -180 and 180 degrees."
            raise ValueError(msg)

        self.latitude = float(latitude)
        self.longitude = float(longitude)
        self.timezone = timezone

        self._check_cache(force_update)

    def _check_cache(self, force_update: bool = False):
        # Check for existence of a cache file
        cache_file = self._find_cache_file(self.latitude, self.longitude)
        if cache_file is None or force_update is True:
            msg = "No cache file or forced update, getting data from OpenMeteo Power."
            self.logger.debug(msg)
            # No cache file, we really have to get the data from the open-meteo server
            self._fetch_data(self.start_date)
            return

        # get age of cache file, if age < 90 days then try to load it. If loading fails retrieve data
        # from the OpenMeteo server .
        r = os.stat(cache_file)
        cache_file_date = datetime.date.fromtimestamp(r.st_mtime)
        age = (datetime.date.today() - cache_file_date).days
        if age < 90:
            msg = "Start loading weather data from cache file: %s" % cache_file
            self.logger.debug(msg)

            status = self._load_cache_file()
            if status is not True:
                msg = "Loading cache file failed, reloading data from OpenMeteo."
                self.logger.debug(msg)
                # Loading cache file failed!
                self._fetch_data(self.start_date)
        else:
            # Cache file is too old. Try loading new data from OpenMeteo
            try:
                msg = "Cache file older then 90 days, reloading data from OpenMeteo."
                self.logger.debug(msg)
                self._fetch_data(self.start_date)
            except Exception as e:
                msg = ("Reloading data from OpenMeteo failed, reverting to (outdated) " +
                       "cache file")
                self.logger.debug(msg)
                status = self._load_cache_file()
                if status is not True:
                    msg = "Outdated cache file failed loading."
                    raise PCSEError(msg)

    def _fetch_data(self, start_date):
        """
        Internal method to fetch and prepare weather data for a given date range.
        Returns a cache file.
        """

        url = self._get_url(previous_runs=True)
        params = {
            "latitude": self.latitude,
            "longitude": self.longitude,
            "start_date": self.format_date(start_date),
            "end_date": self.format_date(self._get_end_date()),
            "daily": self.daily_variables,
            "hourly": self.hourly_variables,
            "timezone": self.timezone,
            "model": self.model,
        }

        response_json = self._get_response(url, params)

        self.elevation = response_json.get('elevation', 0)

        # raw_data = self._extract_weather_data(weather_api_object, params)

        df = self._prepare_weather_dataframe(response_json)

        self._make_WeatherDataContainers(df.to_dict(orient='records'))

        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        self._dump(cache_filename)

    # request routine with error checks
    def _get_response(self, url, params):
        try:
            response = requests.get(url, params=params, timeout=10)
            self._check_response_status(response)

            response_json = response.json()
            return response_json
        except requests.RequestException as e:
            print(f"Error fetching weather data: {e}")
            return None

    def _check_response_status(self, response):
        if response.status_code != self.HTTP_OK:
            msg = ("Failed retrieving OpenMeteo data, server returned HTTP " +
                   "code: %i on following URL %s") % (response.status_code, response.url)
            raise PCSEError(msg)

    def _find_cache_file(self, latitude, longitude):
        """Try to find a cache file for given latitude/longitude.

        Returns None if the cache file does not exist, else it returns the full path
        to the cache file.
        """
        cache_filename = self._get_cache_filename(latitude, longitude)
        if os.path.exists(cache_filename):
            return cache_filename
        else:
            return None

    def _get_cache_filename(self, latitude, longitude):
        """Constructs the filename used for cache files given latitude and longitude

        The latitude and longitude is coded into the filename by truncating on
        0.1 degree. So the cache filename for a point with lat/lon 52.56/-124.78 and using the
        "knmi_seamless" model will be: OpenMeteoWeatherDataProvider_LAT00525_LON-1247_knmi_.cache
        """

        fname = "%s_LAT%05i_LON%05i_%s.cache" % (self.__class__.__name__,
                                                 int(latitude * 10), int(longitude * 10),
                                                 self.model[:5])
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, fname)
        return cache_filename

    def _load_cache_file(self):
        """Loads the data from the cache file. Return True if successful.
        """
        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        try:
            self._load(cache_filename)
            msg = "Cache file successfully loaded."
            self.logger.debug(msg)
            return True
        except (IOError, EnvironmentError, EOFError) as e:
            msg = "Failed to load cache from file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
            return False

    def _make_WeatherDataContainers(self, recs):
        """Create a WeatherDataContainers from recs, compute ET and store the WDC's.
        """

        for rec in recs:
            # Build weather data container from dict 't'
            wdc = WeatherDataContainer(**rec)

            # add wdc to dictionary for this date
            self._store_WeatherDataContainer(wdc, wdc.DAY)

    def _prepare_weather_dataframe(self, weather_data):
        """
        Converts raw Open-Meteo weather data into a single daily DataFrame
        Currently, it is tailored to the inputs required by PCSE

        Returns:
            DataFrame: Daily weather data with dates as index.
        """
        #  Process daily data
        daily = weather_data.get('daily', {})
        df_daily = pd.DataFrame(daily)
        # Convert the 'time' column to datetime objects and set as index
        df_daily['date'] = pd.to_datetime(df_daily['time'])
        df_daily.set_index('date', inplace=True)

        # Rename daily columns for clarity
        df_daily.rename(columns={
            'temperature_2m_min': 'TMIN',
            'temperature_2m_max': 'TMAX',
            'precipitation_sum': 'RAIN',  # in mm/day
            'shortwave_radiation_sum': 'IRRAD'  # in MJ/m²/day
        }, inplace=True)

        #  Process hourly data
        hourly = weather_data.get('hourly', {})
        df_hourly = pd.DataFrame(hourly)
        # Convert hourly time to datetime objects
        df_hourly['date'] = pd.to_datetime(df_hourly['time'])
        # Set time as the DataFrame index
        df_hourly.set_index('date', inplace=True)

        # Compute daily averages from hourly data:
        df_hourly.drop(columns=['time'], inplace=True)
        df_hourly_daily = df_hourly.groupby(df_hourly.index.date).mean()
        # Convert the index back to datetime
        df_hourly_daily.index = pd.to_datetime(df_hourly_daily.index)
        df_hourly_daily.rename(columns={
            'temperature_2m': 'TEMP',
            "wind_speed_10m": 'WIND',
            'dewpoint_2m': 'dewpoint'
        }, inplace=True)

        # Merge on the date index (inner join to keep only days that exist in both)
        df_openmeteo = pd.merge(df_daily, df_hourly_daily, left_index=True, right_index=True, how='inner')

        # Convert irradiation from MJ/m²/day to W/m²/day.
        df_openmeteo['IRRAD'] = df_openmeteo['IRRAD'] * 1e6

        # Convert precipitation from mm/day to cm/day.
        df_openmeteo['RAIN'] = df_openmeteo['RAIN'] * 0.1

        # Convert wind from 10m to 2m
        df_openmeteo['WIND'] = wind10to2(df_openmeteo['WIND'])

        # Calculate vapor pressure (in hPa) from dewpoint (°C) using the formula:
        # e = 6.108 * exp((17.27 * T_d) / (T_d + 237.3))
        df_openmeteo['VAP'] = (6.108 * np.exp((17.27 * df_openmeteo['dewpoint']) / (df_openmeteo['dewpoint'] + 237.3)))

        df_openmeteo.drop(columns=['dewpoint'], inplace=True)

        df_openmeteo['DAY'] = df_openmeteo.index.date

        df_openmeteo['LAT'] = self.latitude
        df_openmeteo['LON'] = self.longitude
        df_openmeteo['ELEV'] = self.elevation

        self.angstA, self.angstB = self._estimate_AngstAB(df_openmeteo)

        df_openmeteo = df_openmeteo[['TMIN', 'TMAX', 'TEMP', 'IRRAD', 'RAIN', 'WIND', 'VAP', 'DAY', 'LAT', 'LON', 'ELEV']]

        E0_list = []
        ES0_list = []
        ET0_list = []

        for row in df_openmeteo.itertuples():
            E0, ES0, ET0 = reference_ET(row.DAY, row.LAT, row.ELEV, row.TMIN,
                                        row.TMAX, row.IRRAD,
                                        row.VAP, row.WIND,
                                        self.angstA, self.angstB, self.ETmodel)

            #  convert to cm/day
            E0_list.append(E0 / 10.)
            ES0_list.append(ES0 / 10.)
            ET0_list.append(ET0 / 10.)

        # Some warning about copying a slice
        df_openmeteo = df_openmeteo.copy()

        df_openmeteo.loc[:, "E0"] = E0_list
        df_openmeteo.loc[:, "ES0"] = ES0_list
        df_openmeteo.loc[:, "ET0"] = ET0_list

        df_openmeteo = df_openmeteo[
            [
                'TMIN',
                'TMAX',
                'TEMP',
                'IRRAD',
                'RAIN',
                'WIND',
                'VAP',
                'DAY',
                'LAT',
                'LON',
                'ELEV',
                'E0',
                'ES0',
                'ET0'
            ]
        ]

        return df_openmeteo

    def calculate_toa_radiation(self, day_of_year):
        """
        Calculate daily Top-of-Atmosphere shortwave radiation
        This ToA estimation was taken from the FAO-56 paper.
        """
        G_sc = 1361  # Solar constant (W/m²)

        # Earth-Sun distance correction factor
        d_r = 1 + 0.033 * np.cos(2 * np.pi * day_of_year / 365)

        # Solar declination (radians)
        delta = np.radians(23.45 * np.sin(2 * np.pi * (day_of_year - 81) / 365))

        # Convert latitude to radians
        phi = np.radians(self.latitude)

        # Sunset hour angle
        h_s = np.arccos(-np.tan(phi) * np.tan(delta))

        # Updated TOA daily radiation (H0)
        H0 = (24 * 3600 / np.pi) * G_sc * d_r * (
                np.cos(phi) * np.cos(delta) * np.sin(h_s) + (h_s * np.sin(phi) * np.sin(delta))
        )

        # Convert from J/m²/day to MJ/m²/day
        H0 = H0 / 1e6

        return H0

    def _estimate_AngstAB(self, df):
        """Determine Angstrom A/B parameters from Top-of-Atmosphere estimation and
        top-of-Canopy (ALLSKY_SFC_SW_DWN) radiation values.

        :param df: dataframe with Openmeteo data
        :return: tuple of Angstrom A/B values

        The Angstrom A/B parameters are determined by dividing swv_dwn by toa_dwn
        and taking the 0.05 percentile for Angstrom A and the 0.98 percentile for
        Angstrom A+B: toa_dwn*(A+B) approaches the upper envelope while
        toa_dwn*A approaches the lower envelope of the records of swv_dwn
        values.

        From PCSE's NASA POWER implementation
        """

        msg = "Start estimation of Angstrom A/B values from Open Meteo data."
        self.logger.debug(msg)

        # check if sufficient data is available to make a reasonable estimate:
        # We want to have at least 200 days available
        if len(df) < 200:
            msg = ("Less then 200 days of data available. Reverting to " +
                   "default Angstrom A/B coefficients (%f, %f)")
            self.logger.warn(msg % (self.angstA, self.angstB))
            return self.angstA, self.angstB

        # calculate relative radiation (swv_dwn/toa_dwn) and percentiles
        doys = pd.to_datetime(df.DAY).dt.dayofyear
        relative_radiation = (df.IRRAD/1e6)/self.calculate_toa_radiation(doys)
        ix = relative_radiation.notnull()
        angstrom_a = float(np.percentile(relative_radiation[ix].values, 5))
        angstrom_ab = float(np.percentile(relative_radiation[ix].values, 98))
        angstrom_b = angstrom_ab - angstrom_a

        try:
            check_angstromAB(angstrom_a, angstrom_b)
        except PCSEError as e:
            msg = ("Angstrom A/B values (%f, %f) outside valid range: %s. " +
                   "Reverting to default values.")
            msg = msg % (angstrom_a, angstrom_b, e)
            self.logger.warn(msg)
            return self.angstA, self.angstB

        msg = "Angstrom A/B values estimated: (%f, %f)." % (angstrom_a, angstrom_b)
        self.logger.debug(msg)

        return angstrom_a, angstrom_b


    def _get_url(self, previous_runs: bool = False) -> str:
        if (self.model in self.dict_forecast_models or self.is_forecast is True) and previous_runs is True:
            return "https://previous-runs-api.open-meteo.com/v1/forecast"
        elif (self.model in self.dict_forecast_models or self.is_forecast is True) and previous_runs is False:
            return "https://api.open-meteo.com/v1/forecast"
        elif self.model in self.dict_historical_models or self.is_forecast is False:
            return "https://archive-api.open-meteo.com/v1/archive"
        else:
            raise ValueError("Model not found. Check model availability.")


    def _get_end_date(self):
        if self.model in self.dict_forecast_models or self.is_forecast is True:
            return datetime.date.today() + datetime.timedelta(days=7)
        elif self.model in self.dict_historical_models or self.is_forecast is False:
            return datetime.date.today() - datetime.timedelta(days=self.delay_historical_models[self.model])
        else:
            raise ValueError("Model not found. Check model availability.")

    def _check_start_date(self):
        if self.start_date is None and self.model in self.dict_forecast_models:
            self.start_date = self.dict_forecast_models[self.model]
        elif self.start_date is None and self.model in self.dict_historical_models:
            self.start_date = self.dict_historical_models[self.model]
        elif self.start_date is None and self.model == "best_match" and self.is_forecast is False:
            self.start_date = self.dict_historical_models["era5"]
        elif self.start_date is None and self.model == "best_match" and self.is_forecast is True:
            self.start_date = self.dict_forecast_models["icon_seamless"]

    @staticmethod
    def format_date(date: Union[str, datetime.date]):
        """
        Converts a date or datetime object to a string in 'YYYY-MM-DD' format.
        If d is already a string, it is returned unchanged.
        """
        if isinstance(date, (datetime.date, datetime)):
            return date.strftime("%Y-%m-%d")
        return date

    @property
    def daily_variables(self):
        return [
        "temperature_2m_max",
        "temperature_2m_min",
        "precipitation_sum",
        "shortwave_radiation_sum",
    ]

    @property
    def hourly_variables(self):
        return [
            "wind_speed_10m",
            "temperature_2m",
            "dewpoint_2m"
        ]


if __name__ == '__main__':
    # Example of grabbing weather from Wageningen
    omwp = OpenMeteoWeatherDataProvider(51.98, 5.65)

    # Get weather for a single day.
    single_date = datetime.date(2024, 5, 15)
    weather_single = omwp(single_date)
    print(f"Weather on {single_date}:", weather_single)
