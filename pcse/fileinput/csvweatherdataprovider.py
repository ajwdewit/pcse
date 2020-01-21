#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl)
# and Zacharias Steinmetz (stei4785@uni-landau.de), Aug 2015
"""A weather data provider reading its data from CSV files.
"""
import os
import datetime as dt
import csv
import math

from ast import literal_eval

from ..base import WeatherDataContainer, WeatherDataProvider
from ..util import reference_ET, angstrom, check_angstromAB
from ..exceptions import PCSEError
from ..settings import settings


class ParseError(PCSEError):
    pass


class OutOfRange(PCSEError):
    pass


class IRRADFromSunshineDuration:

    def __init__(self, latitude, angstA, angstB):

        assert -90 < latitude < 90, \
            "Invalid latitude value (%s) encountered" % latitude
        check_angstromAB(angstA, angstB)
        self.latitude = latitude
        self.angstA = angstA
        self.angstB = angstB

    def __call__(self, value, day):
        """Computes irradiance in J/m2/day from sunshine duration by applying the Angstrom equation

        :param value: sunshine duration in hours
        :param day: the day
        :return: irradiance in J/m2/day
        """
        assert 0 <= value <= 24, \
            "Invalid sunshine duration value (%s) encountered at day %s" % (value, day)
        irrad = angstrom(day, self.latitude, value, self.angstA, self.angstB)

        return irrad


def csvdate_to_date(x, dateformat):
    """Converts string x to a datetime.date using given format.

    :param x: the string representing a date
    :param dateformat: a strptime() accepted date format
    :return: a date
    """
    return dt.datetime.strptime(x, dateformat).date()


# Conversion functions
def NoConversion(x, d):
    return float(x)


def kJ_to_J(x, d):
    return float(x)*1000.


def mm_to_cm(x, d):
    return float(x)/10.


def kPa_to_hPa(x, d):
    return float(x)*10.


class CSVWeatherDataProvider(WeatherDataProvider):
    """Reading weather data from a CSV file.

    :param csv_fname: name of the CSV file to be read
    :param delimiter: CSV delimiter
    :param dateformat: date format to be read. Default is '%Y%m%d'
    :keyword ETmodel: "PM"|"P" for selecting Penman-Monteith or Penman
        method for reference evapotranspiration. Default is 'PM'.
    :param force_reload: Ignore cache file and reload from the CSV file

    The CSV file should have the following structure (sample), missing values should be added as 'NaN'::

        ## Site Characteristics
        Country     = 'Netherlands'
        Station     = 'Wageningen, Haarweg'
        Description = 'Observed data from Station Haarweg in Wageningen'
        Source      = 'Meteorology and Air Quality Group, Wageningen University'
        Contact     = 'Peter Uithol'
        Longitude = 5.67; Latitude = 51.97; Elevation = 7; AngstromA = 0.18; AngstromB = 0.55; HasSunshine = False
        ## Daily weather observations (missing values are NaN)
        DAY,IRRAD,TMIN,TMAX,VAP,WIND,RAIN,SNOWDEPTH
        20040101,NaN,-0.7,1.1,0.55,3.6,0.5,NaN
        20040102,3888,-7.5,0.9,0.44,3.1,0,NaN
        20040103,2074,-6.8,-0.5,0.45,1.8,0,NaN
        20040104,1814,-3.6,5.9,0.66,3.2,2.5,NaN
        20040105,1469,3,5.7,0.78,2.3,1.3,NaN
        [...]

        with
        IRRAD in kJ/m2/day or hours
        TMIN and TMAX in Celsius (°C)
        VAP in kPa
        WIND in m/sec
        RAIN in mm
        SNOWDEPTH in cm

    For reading weather data from a file, initially the CABOWeatherDataProvider
    was available which read its data from text in the CABO weather format.
    Nevertheless, building CABO weather files is tedious as for each year a new
    file must constructed. Moreover it is rather error prone and formatting
    mistakes are easily leading to errors.

    To simplify providing weather data to PCSE models, a new data provider
    has been derived from the ExcelWeatherDataProvider that reads its data
    from simple CSV files.

    The CSVWeatherDataProvider assumes that records are complete and does
    not make an effort to interpolate data as this can be easily
    accomplished in a text editor. Only SNOWDEPTH is allowed to be missing
    as this parameter is usually not provided outside the winter season.
    """

    obs_conversions = {
        "TMAX": NoConversion,
        "TMIN": NoConversion,
        "IRRAD": kJ_to_J,
        "VAP": kPa_to_hPa,
        "WIND": NoConversion,
        "RAIN": mm_to_cm,
        "SNOWDEPTH": NoConversion
    }

    def __init__(self, csv_fname, delimiter=',', dateformat='%Y%m%d',
                 ETmodel='PM', force_reload=False):
        WeatherDataProvider.__init__(self)

        self.fp_csv_fname = os.path.abspath(csv_fname)
        self.dateformat = dateformat
        self.ETmodel = ETmodel
        if not os.path.exists(self.fp_csv_fname):
            msg = "Cannot find weather file at: %s" % self.fp_csv_fname
            raise PCSEError(msg)

        if force_reload or not self._load_cache_file(self.fp_csv_fname):
            with open(csv_fname, 'r') as csv_file:
                self._read_meta(csv_file)
                self._read_observations(csv_file, delimiter)
            self._write_cache_file(self.fp_csv_fname)

    def _read_meta(self, csv_file):
        header = {}
        for line in csv_file:
            if line.startswith('## Daily weather observations'):
                break
            statements = line.split(';')
            for stmt in statements:
                key, val = stmt.split('=')
                header[key.strip()] = literal_eval(val.strip())

        self.nodata_value = -99
        self.description = [u"Weather data for:",
                            u"Country: %s" % header['Country'],
                            u"Station: %s" % header['Station'],
                            u"Description: %s" % header['Description'],
                            u"Source: %s" % header['Source'],
                            u"Contact: %s" % header['Contact']]

        self.longitude = float(header['Longitude'])
        self.latitude = float(header['Latitude'])
        self.elevation = float(header['Elevation'])
        angstA = float(header['AngstromA'])
        angstB = float(header['AngstromB'])
        self.angstA, self.angstB = check_angstromAB(angstA, angstB)
        self.has_sunshine = bool(header['HasSunshine'])

        # If the file has sunshine duration, we replace the convertor with the angstrom module
        if self.has_sunshine:
            self.obs_conversions["IRRAD"] = IRRADFromSunshineDuration(self.latitude, self.angstA, self.angstB)

    def _read_observations(self, csv_file, delimiter):
        """Processes the rows with meteo data and converts into the correct units.
        """
        obs = csv.DictReader(csv_file, delimiter=delimiter, quotechar='"')
        for i, d in enumerate(obs):
            try:
                day = None
                day = csvdate_to_date(d.pop("DAY"), self.dateformat)
                row = {"DAY":  day}
                for label in self.obs_conversions.keys():
                    func = self.obs_conversions[label]
                    value = float(d[label])
                    r = func(value, day)
                    if math.isnan(r):
                        if label == "SNOWDEPTH":
                            continue
                        raise ParseError
                    row[label] = r

                # Reference ET in mm/day
                e0, es0, et0 = reference_ET(LAT=self.latitude, ELEV=self.elevation,
                                            ANGSTA=self.angstA, ANGSTB=self.angstB,
                                            ETMODEL=self.ETmodel, **row)
                # convert to cm/day
                row["E0"] = e0/10.
                row["ES0"] = es0/10.
                row["ET0"] = et0/10.

                wdc = WeatherDataContainer(LAT=self.latitude, LON=self.longitude, ELEV=self.elevation, **row)
                self._store_WeatherDataContainer(wdc, day)
            except (ParseError, KeyError):
                msg = "Failed reading element '%s' for day '%s' at line %i. Skipping ..." % (label, day, i)
                self.logger.warn(msg)
            except ValueError as e:  # strange value in cell
                msg = "Failed computing a value for day '%s' at row %i" % (day, i)
                self.logger.warn(msg)

    def _load_cache_file(self, csv_fname):

        cache_filename = self._find_cache_file(csv_fname)
        if cache_filename is None:
            return False
        else:
            self._load(cache_filename)
            return True

    def _find_cache_file(self, csv_fname):
        """Try to find a cache file for file name

        Returns None if the cache file does not exist, else it returns the full
        path to the cache file.
        """
        cache_filename = self._get_cache_filename(csv_fname)
        if os.path.exists(cache_filename):
            cache_date = os.stat(cache_filename).st_mtime
            csv_date = os.stat(csv_fname).st_mtime
            if cache_date > csv_date:  # cache is more recent then CSV file
                return cache_filename

        return None

    def _get_cache_filename(self, csv_fname):
        """Constructs the filename used for cache files given csv_fname
        """
        basename = os.path.basename(csv_fname)
        filename, ext = os.path.splitext(basename)

        tmp = "%s_%s.cache" % (self.__class__.__name__, filename)
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, tmp)
        return cache_filename

    def _write_cache_file(self, csv_fname):

        cache_filename = self._get_cache_filename(csv_fname)
        try:
            self._dump(cache_filename)
        except (IOError, EnvironmentError) as e:
            msg = "Failed to write cache to file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
