# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import os, sys
import glob
import calendar
import numpy as np
import datetime as dt
import warnings
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

def xlsdate_to_date(value, sheet):
    """Convert an excel date into a python date

    :param value: A value from an excel cell
    :return: a python date
    """
    year, month, day, hr, min, sec = xlrd.xldate_as_tuple(value, sheet.book.datemode)
    return dt.date(year, month, day)

class ExcelWeatherDataProvider(WeatherDataProvider):
    """Reading weather data from Excel (xls & xlst) files.
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
    site_row = 8
    label_row = 10
    data_start_row = 12

    def __init__(self, fname):
        WeatherDataProvider.__init__(self)

        fp_fname = os.path.abspath(fname)
        if not os.path.exists(fname):
            msg = "Cannot find weather file at: %s" % fp_fname
            raise PCSEError(msg)

        book = xlrd.open_workbook(fp_fname)
        sheet = book.sheet_by_index(0)

        self._read_header(sheet)
        self._read_site_characteristics(sheet)
        self._read_observations(sheet)

    def _read_header(self, sheet):

        country = sheet.cell_value(1, 1)
        station = sheet.cell_value(2, 1)
        desc = sheet.cell_value(3, 1)
        src = sheet.cell_value(4, 1)
        contact = sheet.cell_value(5, 1)
        self.nodata_value = sheet.cell_value(6, 1)
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

        labels = [cell.value for cell in sheet.row(self.label_row)]
        # Start reading all rows with data
        rownums = range(sheet.nrows)
        for rownum in rownums[self.data_start_row:]:
            row = sheet.row(rownum)
            d = {}
            nextrow = False
            for cell, label in zip(row, labels):
                try:
                    value = float(cell.value)
                except ValueError as e:
                    msg = "Failed reading row: %i. Skipping..." % (rownum + 1)
                    self.logger.warning(msg)
                    print(msg)
                    nextrow = True

                if label == "IRRAD" and self.has_sunshine is True:
                    d[label] = angstrom(d["DAY"], self.latitude, float(cell.value),
                                        self.angstA, self.angstB)
                func = self.obs_conversions[label]
                d[label] = func(float(cell.value), sheet)
            e0, es0, et0 = reference_ET(LAT=self.latitude, ELEV=self.elevation, ANGSTA=self.angstA,
                                        ANGSTB=self.angstB, **d)
            d["E0"] = e0; d["ES0"] = es0; d["ET0"] = et0
            wdc = WeatherDataContainer(LAT=self.latitude, LON=self.longitude, ELEV=self.elevation, **d)
            self._store_WeatherDataContainer(wdc, d["DAY"])

