# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2015
"""A weather data provider reading its data from Excel files.
"""
import os
import datetime as dt
import xlrd

from ..base_classes import WeatherDataContainer, WeatherDataProvider
from ..util import reference_ET, angstrom, check_angstromAB
from ..exceptions import PCSEError
from ..settings import settings

# Conversion functions, note that we defined a second parameter 's' that is only there
# to catch the sheet which is needed in the function xlsdate_to_date(value, sheet).
NoConversion = lambda x, s: x
kJ_to_J = lambda x, s: x*1000.
kPa_to_hPa = lambda x, s: x*10.
mm_to_cm = lambda x, s: x/10.


class NoDataError(PCSEError):
    pass


class OutOfRange(PCSEError):
    pass


def xlsdate_to_date(value, sheet):
    """Convert an excel date into a python date

    :param value: A value from an excel cell
    :param sheet: A reference to the excel sheet for getting the datemode
    :return: a python date
    """
    year, month, day, hr, min, sec = xlrd.xldate_as_tuple(value, sheet.book.datemode)
    return dt.date(year, month, day)


class ExcelWeatherDataProvider(WeatherDataProvider):
    """Reading weather data from an excel file.

    :param xls_fname: name of the Excel file to be read
    :param mising_snow_depth: the value that should use for missing SNOW_DEPTH values

    For reading weather data from file, initially only the CABOWeatherDataProvider
    was available that reads its data from a text file in the CABO Weather format.
    Nevertheless, building CABO weather files is tedious as for each year a new
    file must constructed. Moreover it is rather error prone and formatting
    mistakes are easily leading to errors.

    To simplify providing weather data to PCSE models, a new data provider
    was written that reads its data from simple excel files

    The ExcelWeatherDataProvider assumes that records are complete and does
    not make an effort to interpolate data as this can be easily
    accomplished in Excel itself. Only SNOW_DEPTH is allowed to be missing
    as this parameter is usually not provided outside the winter season.
    """
    obs_conversions = {
        "TMAX": NoConversion,
        "TMIN": NoConversion,
        "IRRAD": kJ_to_J,
        "DAY": xlsdate_to_date,
        "VAP": kPa_to_hPa,
        "WIND": NoConversion,
        "RAIN": mm_to_cm,
        "SNOWDEPTH": NoConversion
    }

    # row numbers where values start. Note that the row numbers are
    # zero-based, so add 1 to find the corresponding row in excel.
    site_row = 8
    label_row = 10
    data_start_row = 12

    def __init__(self, xls_fname, missing_snow_depth=None):
        WeatherDataProvider.__init__(self)

        self.fp_xls_fname = os.path.abspath(xls_fname)
        self.missing_snow_depth = missing_snow_depth
        if not os.path.exists(self.fp_xls_fname):
            msg = "Cannot find weather file at: %s" % self.fp_xls_fname
            raise PCSEError(msg)

        if not self._load_cache_file(self.fp_xls_fname):  # Cache file cannot be loaded
            book = xlrd.open_workbook(self.fp_xls_fname)
            sheet = book.sheet_by_index(0)

            self._read_header(sheet)
            self._read_site_characteristics(sheet)
            self._read_observations(sheet)

            self._write_cache_file(self.fp_xls_fname)

    def _read_header(self, sheet):

        country = sheet.cell_value(1, 1)
        station = sheet.cell_value(2, 1)
        desc = sheet.cell_value(3, 1)
        src = sheet.cell_value(4, 1)
        contact = sheet.cell_value(5, 1)
        self.nodata_value = float(sheet.cell_value(6, 1))
        self.description = [u"Weather data for:",
                            u"Country: %s" % country,
                            u"Station: %s" % station,
                            u"Description: %s" % desc,
                            u"Source: %s" % src,
                            u"Contact: %s" % contact]

    def _read_site_characteristics(self, sheet):

        self.longitude = float(sheet.cell_value(self.site_row, 0))
        self.latitude = float(sheet.cell_value(self.site_row, 1))
        self.elevation = float(sheet.cell_value(self.site_row, 2))
        angstA = float(sheet.cell_value(self.site_row, 3))
        angstB = float(sheet.cell_value(self.site_row, 4))
        self.angstA, self.angstB = check_angstromAB(angstA, angstB)
        self.has_sunshine = bool(sheet.cell_value(self.site_row, 5))

    def _read_observations(self, sheet):

        # First get the column labels
        labels = [cell.value for cell in sheet.row(self.label_row)]

        # Start reading all rows with data
        rownums = list(range(sheet.nrows))
        for rownum in rownums[self.data_start_row:]:
            try:
                row = sheet.row(rownum)
                d = {}
                for cell, label in zip(row, labels):
                    # explicitly convert to float. If this fails a ValueError will be thrown
                    value = float(cell.value)

                    # Check for observations marked as missing. Currently only missing
                    # data is allowed for SNOWDEPTH. Otherwise raise an error
                    eps = 0.0001
                    if (value - self.nodata_value) < eps:
                        if label == "SNOWDEPTH":
                            value = self.missing_snow_depth
                        else:
                            raise NoDataError

                    if label == "IRRAD" and self.has_sunshine is True:
                        if 0 < value < 24:
                            d[label] = angstrom(d["DAY"], self.latitude, value, self.angstA, self.angstB)
                        else:
                            msg = "Sunshine duration not within 0-24 interval for row %i" % (rownum + 1)
                            raise OutOfRange(msg)

                    func = self.obs_conversions[label]
                    d[label] = func(value, sheet)

                # Reference ET in mm/day
                e0, es0, et0 = reference_ET(LAT=self.latitude, ELEV=self.elevation, ANGSTA=self.angstA,
                                            ANGSTB=self.angstB, **d)
                # convert to cm/day
                d["E0"] = e0/10.; d["ES0"] = es0/10.; d["ET0"] = et0/10.

                wdc = WeatherDataContainer(LAT=self.latitude, LON=self.longitude, ELEV=self.elevation, **d)
                self._store_WeatherDataContainer(wdc, d["DAY"])

            except ValueError as e: # strange value in cell
                msg = "Failed reading row: %i. Skipping..." % (rownum + 1)
                self.logger.warn(msg)
                print(msg)

            except NoDataError as e: # Missing value encountered
                msg = "No data value (%f) encountered at row %i. Skipping..." % (self.nodata_value, (rownum + 1))
                self.logger.warn(msg)

            except OutOfRange as e:
                self.logger.warn(e.message)

    def _load_cache_file(self, xls_fname):

        cache_filename = self._find_cache_file(xls_fname)
        if cache_filename is None:
            return False
        else:
            self._load(cache_filename)
            return True

    def _find_cache_file(self, xls_fname):
        """Try to find a cache file for file name

        Returns None if the cache file does not exist, else it returns the full path
        to the cache file.
        """
        cache_filename = self._get_cache_filename(xls_fname)
        if os.path.exists(cache_filename):
            cache_date = os.stat(cache_filename).st_mtime
            xls_date = os.stat(xls_fname).st_mtime
            if cache_date > xls_date:  # cache is more recent then XLS file
                return cache_filename

        return None

    def _get_cache_filename(self, xls_fname):
        """Constructs the filename used for cache files given xls_fname
        """
        basename = os.path.basename(xls_fname)
        filename, ext = os.path.splitext(basename)

        tmp = "%s_%s.cache" % (self.__class__.__name__, filename)
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, tmp)
        return cache_filename

    def _write_cache_file(self, xls_fname):

        cache_filename = self._get_cache_filename(xls_fname)
        try:
            self._dump(cache_filename)
        except (IOError, EnvironmentError) as e:
            msg = "Failed to write cache to file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
