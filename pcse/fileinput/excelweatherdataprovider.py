# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2015
"""A weather data provider reading its data from Excel files.
"""
import os
import openpyxl

from ..base import WeatherDataContainer, WeatherDataProvider
from ..util import reference_ET, angstrom, check_angstromAB
from ..exceptions import PCSEError
from ..settings import settings

# Conversion functions
NoConversion = lambda x: x
kJ_to_J = lambda x: x*1000.
kPa_to_hPa = lambda x: x*10.
mm_to_cm = lambda x: x/10.


def determine_true_false(value):
    """OpenPyXL has a somewhat strange treatment of true/false

    Excel cell value      OpenPyXL cell value       Type
       =FALSE()            '=FALSE()'               str
       =FALSE              '=FALSE'                 str
       false                False                   bool
       =TRUE()             '=TRUE()'                str
       =TRUE               '=TRUE'                  str
       TRUE                 True                    bool

    Finally, there is the possibility for true/false to be presented
    by 0/1 integer values.

    This function tries to handle all of the them.
    """
    if isinstance(value, str):
        v = value.lower()
        if "true" in v:
            return True
        elif "false" in v:
            return False
    elif isinstance(value, bool):
        return value
    elif isinstance(value, int):
        if value == 1:
            return True
        elif value == 0:
            return False

    msg = f"cannot determine True|False: {value}"
    raise ValueError(msg)


class ExcelWeatherDataProvider(WeatherDataProvider):
    """Reading weather data from an excel file (.xlsx only).

    :param xls_fname: name of the Excel file to be read
    :param mising_snow_depth: the value that should use for missing SNOW_DEPTH values,
           the default value is `None`.
    :param force_reload: bypass the cache file, reload data from the .xlsx file and
           write a new cache file. Cache files are written under `$HOME/.pcse/meteo_cache`

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
        "VAP": kPa_to_hPa,
        "WIND": NoConversion,
        "RAIN": mm_to_cm,
        "SNOWDEPTH": NoConversion
    }

    # row numbers where values start. Note that the row numbers are
    # zero-based, so add 1 to find the corresponding row in excel.
    site_row = 9
    label_row = 11
    data_start_row = 13

    def __init__(self, xls_fname, missing_snow_depth=None, force_reload=False):
        WeatherDataProvider.__init__(self)

        self.fp_xls_fname = os.path.abspath(xls_fname)
        self.missing_snow_depth = missing_snow_depth
        if not os.path.exists(self.fp_xls_fname):
            msg = "Cannot find weather file at: %s" % self.fp_xls_fname
            raise PCSEError(msg)

        if force_reload or not self._load_cache_file(self.fp_xls_fname):
            book = openpyxl.load_workbook(self.fp_xls_fname, read_only=True)
            sheet = book.active

            self._read_header(sheet)
            self._read_site_characteristics(sheet)
            self._read_observations(sheet)

            self._write_cache_file(self.fp_xls_fname)

    def _read_header(self, sheet):

        country = sheet["B2"].value
        station = sheet["B3"].value
        desc = sheet["B4"].value
        src = sheet["B5"].value
        contact = sheet["B6"].value
        self.nodata_value = float(sheet["B7"].value)
        self.description = [u"Weather data for:",
                            u"Country: %s" % country,
                            u"Station: %s" % station,
                            u"Description: %s" % desc,
                            u"Source: %s" % src,
                            u"Contact: %s" % contact]

    def _read_site_characteristics(self, sheet):

        self.longitude = float(sheet[f"A{self.site_row}"].value)
        self.latitude = float(sheet[f"B{self.site_row}"].value)
        self.elevation = float(sheet[f"C{self.site_row}"].value)
        angstA = float(sheet[f"D{self.site_row}"].value)
        angstB = float(sheet[f"E{self.site_row}"].value)
        self.angstA, self.angstB = check_angstromAB(angstA, angstB)
        try:
            has_sunshine = sheet[f"F{self.site_row}"].value
            self.has_sunshine = determine_true_false(has_sunshine)
        except ValueError as e:
            raise PCSEError(f"Cannot determine if sheet as radiation or sunshine hours: {e}")


    def _read_observations(self, sheet):

        # First get the column labels
        labels = [cell.value for cell in sheet[self.label_row]]

        # Start reading all rows with data
        # rownums = list(range(sheet.nrows))
        for rownum, row in enumerate(sheet[self.data_start_row:sheet.max_row]):
            try:
                d = {}
                for cell, label in zip(row, labels):
                    if label == "DAY":
                        if cell.value is None:
                            raise ValueError
                        else:
                            d[label] = cell.value.date()
                            continue

                    # explicitly convert to float. If this fails a ValueError will be thrown
                    value = float(cell.value)

                    # Check for observations marked as missing. Currently only missing
                    # data is allowed for SNOWDEPTH. Otherwise raise an error
                    if self._is_missing_value(value):
                        if label == "SNOWDEPTH":
                            value = self.missing_snow_depth
                        else:
                            raise ValueError()

                    if label == "IRRAD" and self.has_sunshine is True:
                        if 0 < value < 24:
                            # Use Angstrom equation to convert sunshine duration to radiation in J/m2/day
                            value = angstrom(d["DAY"], self.latitude, value, self.angstA, self.angstB)
                            value /= 1000.  # convert to kJ/m2/day for compatibility with obs_conversion function
                        else:
                            msg = "Sunshine duration not within 0-24 interval for row %i" % \
                                  (rownum + self.data_start_row)
                            raise ValueError(msg)

                    func = self.obs_conversions[label]
                    d[label] = func(value)

                # Reference ET in mm/day
                e0, es0, et0 = reference_ET(LAT=self.latitude, ELEV=self.elevation, ANGSTA=self.angstA,
                                            ANGSTB=self.angstB, **d)
                # convert to cm/day
                d["E0"] = e0/10.; d["ES0"] = es0/10.; d["ET0"] = et0/10.

                wdc = WeatherDataContainer(LAT=self.latitude, LON=self.longitude, ELEV=self.elevation, **d)
                self._store_WeatherDataContainer(wdc, d["DAY"])

            except ValueError as e:  # strange value in cell
                msg = "Failed reading row: %i. Skipping..." % (rownum + self.data_start_row)
                self.logger.warning(msg)
                print(msg)

    def _load_cache_file(self, xls_fname):

         cache_filename = self._find_cache_file(xls_fname)
         if cache_filename is None:
             return False
         else:
             try:
                 self._load(cache_filename)
                 return True
             except:
                 return False

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

    def _is_missing_value(self, value):
        """Checks if value is equal to the value specified for missign date

        :return: True|False
        """
        eps = 0.0001
        if abs(value - self.nodata_value) < eps:
            return True
        else:
            return False