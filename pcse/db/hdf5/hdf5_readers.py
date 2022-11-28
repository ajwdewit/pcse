# -*- coding: utf-8 -*-
# Copyright (c) 2004-2017 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), November 2017
"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS12 database schema but assume
that tables are stored in a pandas HDFstore. See also:
https://pandas.pydata.org/pandas-docs/stable/io.html#io-hdf5
"""
import os, sys
import datetime as dt

import pandas as pd
import numpy as np
import yaml

from ...base import WeatherDataProvider
from ... import exceptions as exc
from ...util import check_date
from ..wofost_parameters import WOFOST_optional_parameters, WOFOST_parameter_codes_single, \
    WOFOST_parameter_codes_tabular, WOFOST_parameters_additional, \
    WOFOST_single2tabular


def fetch_crop_name(HDFstore, crop_no):
    """Retrieves the name of the crop from the CROP table for
    given crop_no.

    :param HDFstore: path to HDFstore to retrieve data from
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """
    HDFstore = os.path.abspath(HDFstore)
    if not os.path.exists(HDFstore):
        msg = "Cannot find HDFstore at: %s" % HDFstore
        raise exc.PCSEError(msg)

    # Read crop name
    wh = "crop_no=%i" % crop_no
    df_crop = pd.read_hdf(HDFstore, "/crop", where=wh)
    if len(df_crop) == 0:
        msg = "Failure reading crop type: crop_no==%i not found!" % crop_no
        raise exc.WeatherDataProviderError(msg)

    return df_crop.crop_name.iloc[0]


class WeatherObsGridDataProvider(WeatherDataProvider):

    def __init__(self, HDFstore=None, grid_no=None):
        WeatherDataProvider.__init__(self)

        if HDFstore is None or grid_no is None:
            msg = "Provide location of HDF weather data store and grid ID."
            raise RuntimeError(msg)

        self.HDFstore = os.path.abspath(HDFstore)
        if not os.path.exists(self.HDFstore):
            msg = "Cannot find HDF5 file at: %s" % self.HDFstore
            raise RuntimeError(msg)

        # Read grid information and meteo time-series
        self.df_grid = pd.read_hdf(self.HDFstore, '/grid', where="grid_no==%i" % grid_no)
        self.df_meteo = pd.read_hdf(self.HDFstore, ('/%s/data' % grid_no))

        if len(self.df_grid) == 0:
            msg = "Cannot find grid_no == %i in grid table." % grid_no
            raise exc.WeatherDataProviderError(msg)

        # Post-processing on meteo records
        self.df_meteo["DAY"] = self.df_meteo.index.date
        self.df_meteo["LAT"] = np.float32(self.df_grid.latitude.iloc[0])
        self.df_meteo["LON"] = np.float32(self.df_grid.longitude.iloc[0])
        self.df_meteo["ELEV"] = np.float32(self.df_grid.altitude.iloc[0])
        self.df_meteo["DTEMP"] = (0.5 * (self.df_meteo.TEMP + self.df_meteo.TMAX)).astype(np.float32)
        self.df_meteo["TMINRA"] = self.df_meteo.TMIN.rolling(7, min_periods=1).mean().astype(np.float32)

        # Metadata for weatherdataprovider
        self.description = "Weather data for grid %i from HDF5 store: %s" % (grid_no, self.HDFstore)
        self.latitude = float(self.df_grid.latitude)
        self.longitude = float(self.df_grid.longitude)
        self.elevation = float(self.df_grid.altitude)
        self._first_date = min(self.df_meteo.DAY)
        self._last_date = max(self.df_meteo.DAY)
        self.store = self.df_meteo.DAY

    def __call__(self, day):

        kday = self.check_keydate(day)
        df = self.df_meteo[self.df_meteo.index == kday.isoformat()]
        if len(df) == 0:
            msg = "Cannot find weather data for %s" % kday
            raise exc.WeatherDataProviderError(msg)
        return list(df.itertuples())[0]

    def export(self):
        return self.df_meteo.to_dict(orient='records')


class CropDataProvider(dict):

    def __init__(self, HDFstore, grid_no, crop_no, campaign_year):

        self.HDFstore = os.path.abspath(HDFstore)
        self.crop_no = int(crop_no)
        self.grid_no = int(grid_no)
        self.campaign_year = int(campaign_year)

        if not os.path.exists(self.HDFstore):
            msg = "Cannot find HDFstore at: %s" % self.HDFstore
            raise exc.WeatherDataProviderError(msg)

        # Read crop name
        self.crop_name = fetch_crop_name(self.HDFstore, self.crop_no)

        # Read crop calendar for variety number
        wh = "grid_no=%i and crop_no=%i and year=%i" % (self.grid_no, self.crop_no, self.campaign_year)
        df_cc = pd.read_hdf(self.HDFstore, "/crop_calendar", where=wh)
        if len(df_cc) == 0:
            msg = "Cannot find crop calendar for: %s" % wh
            raise exc.WeatherDataProviderError(msg)
        self.variety_no = df_cc.variety_no.iloc[0]

        # Read parameter values from crop and variety parameter tables
        self._read_crop_parameters()
        self._read_variety_parameters()
        # Some additional parameters not in CGMS
        self.update(WOFOST_parameters_additional)

        # Finally add crop name
        self["CRPNAM"] = self.crop_name

    def _read_crop_parameters(self):

        wh = "crop_no=%i" % self.crop_no
        df_croppar = pd.read_hdf(self.HDFstore, "/crop_parameter_value", where=wh)

        for parcode in WOFOST_parameter_codes_single:
            df = df_croppar[df_croppar.parameter_code == parcode]
            if len(df) == 0:
                if parcode not in WOFOST_optional_parameters:
                    msg = "Failure retrieving parameter '%s'" % parcode
                    raise exc.WeatherDataProviderError(msg)
            self[parcode] = df.parameter_xvalue.iloc[0]

            if parcode not in WOFOST_single2tabular:
                self[parcode] = df.parameter_xvalue.iloc[0]
            else:
                pvalue = df.parameter_xvalue.iloc[0]
                code, value = self._convert_single2tabular(parcode, pvalue)
                self[code] = value

        for parcode in WOFOST_parameter_codes_tabular:
            df = df_croppar[df_croppar.parameter_code.str.startswith(parcode)]
            if len(df) == 0:
                if parcode not in WOFOST_optional_parameters:
                    msg = "Failure retrieving parameter '%s'" % parcode
                    raise exc.WeatherDataProviderError(msg)
            if len(df) == 1:
                msg = ("Single parameter value found for crop_no=%s, parameter_code='%s' while "
                       "tabular parameter expected." % (self.crop_no, parcode))
                raise exc.WeatherDataProviderError(msg)
            table = []
            for t in df.itertuples():
                table.extend([t.parameter_xvalue, t.parameter_yvalue])
            self[parcode] = table

    def _read_variety_parameters(self):

        wh = "crop_no=%i and variety_no=%i" % (self.crop_no, self.variety_no)
        df_croppar = pd.read_hdf(self.HDFstore, "/variety_parameter_value", where=wh)

        for parcode in WOFOST_parameter_codes_single:
            df = df_croppar[df_croppar.parameter_code == parcode]
            if len(df) == 0:
                continue
            if parcode not in WOFOST_single2tabular:
                self[parcode] = df.parameter_xvalue.iloc[0]
            else:
                pvalue = df.parameter_xvalue.iloc[0]
                code, value = self._convert_single2tabular(parcode, pvalue)
                self[code] = value

        for parcode in WOFOST_parameter_codes_tabular:
            df = df_croppar[df_croppar.parameter_code.str.startswith(parcode)]
            if len(df) == 0:
                continue
            if len(df) == 1:
                msg = ("Single parameter value found for crop_no=%s, parameter_code='%s' while "
                       "tabular parameter expected." % (self.crop_no, parcode))
                raise exc.WeatherDataProviderError(msg)
            table = []
            for t in df.itertuples():
                table.extend([t.parameter_xvalue, t.parameter_yvalue])
            self[parcode] = table

    def _convert_single2tabular(self, parameter_code, pvalue):
        """Converts the single parameter into a tabular parameter.
        """
        tabular_parameter_code, template = WOFOST_single2tabular[parameter_code]
        tabular_values = [pvalue if v is None else v for v in template]

        return tabular_parameter_code, tabular_values


class SoilDataProviderSingleLayer(dict):
    """Class for providing soil data from the ROOTING_DEPTH AND
    SOIL_PHYSICAL_GROUP tableS in a CGMS8/12 database. This
    applies to the single layered soil only.

    :param HDFstore: Path to HDF5 store where tables and data can be found
    :param stu_no: Integer stu no, maps to the STU_NO column in the table
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

    def __init__(self, HDFstore, stu_no):
        dict.__init__(self)

        self.HDFstore = os.path.abspath(HDFstore)

        # First get the rooting depth class and soil_group_no
        rd_class, spg_no = self._get_from_STU(stu_no)
        # Get the actual rooting depth [cm]
        self._get_rooting_depth(rd_class)
        # Get the actual soil hydrological parameters.
        self._get_soil_hydraulic_parameters(spg_no)

        # define SMLIM
        self["SMLIM"] = self["SMFCF"]

    def _get_from_STU(self, stu_no):
        """Retrieves the soil parameters for the given soil typologic unit
        (stu_no) from the tables SOIL_PHYSICAL_GROUP and ROOTING_DEPTH.
        """
        wh = "stu_no = %i" % stu_no
        df = pd.read_hdf(self.HDFstore, "/soil_typologic_unit", where=wh)
        if len(df) == 0:
            msg = ("No record found for stu_no=%i in table SOIL_TYPOLOGIC_UNIT." % stu_no)
            raise exc.PCSEError(msg)
        soil_group_no = int(df.soil_group_no.iloc[0])
        rd_class = int(df.calculated_rooting_depth.iloc[0])

        return (rd_class, soil_group_no)

    def _get_rooting_depth(self, rd_class):
        """Gets the rooting depth from the table ROOTING_DEPTH and
         stores into self[] directly under parameter name 'RDMSOL'.

        :param rd_class: The rooting depth class (integer)
        """

        df = pd.read_hdf(self.HDFstore, "/rooting_depth")
        # Note that we need to loop over the rows of df here instead of putting a
        # where="class = rd_class" in the read_hdf() call because 'class' is a reserved
        # word and we get a SyntaxError if we do otherwise.
        for t in df.itertuples(index=False):
            if t[0] == rd_class:
                break
        else:
            msg = ("No record found for rooting depth class %i in table ROOTING_DEPTH." % rd_class)
            raise exc.PCSEError(msg)

        self["RDMSOL"] = t.min_depth

    def _get_soil_hydraulic_parameters(self, spg_no):
        """Retrieves the soil hydraulic parameters and stores into self[] directly.

        :param spg_no: the soil physical group number (integer)
        """

        wh = "soil_group_no = %i" % spg_no
        df = pd.read_hdf(self.HDFstore, "/soil_physical_group", where=wh)
        for (wofost_soil_par, db_soil_par) in self.soil_parameters:
            dfs = df[df.parameter_code == db_soil_par]
            if len(dfs) == 0:
                msg = "Parameter %s not found in table SOIL_PHYSICAL_GROUP" % db_soil_par
                raise exc.PCSEError(msg)
            self[wofost_soil_par] = dfs.parameter_xvalue.iloc[0]


class AgroManagementDataProvider(list):
    """Class for providing agromanagement data from the CROP_CALENDAR table in a HDF5 file.

    :param HDFstore: path to HDFstore to retrieve data from
    :param grid_no: Integer grid ID, maps to the grid_no column in the table
    :param crop_no: Integer id of crop, maps to the crop_no column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
           The campaign year usually refers to the year of the harvest. Thus for crops
           crossing calendar years, the start_date can be in the previous year.
    :keyword campaign_start: Optional keyword that can be used to define the start of the
           campaign. Note that by default the campaign_start_date is set equal to the
           crop_start_date which means that the simulation starts when the crop starts.
           This default behaviour can be changed using this keyword. It can have multiple meanings:

               - if a date object is passed, the campaign is assumed to start on this date.
               - if an int/float is passed, the campaign_start_date is calculated as the
                 crop_start_date minus the number of days provided by campaign_start.

    For adjusting the campaign_start_Date, see also the `set_campaign_start_date(date)` method
    to update the campaign_start_date on an existing AgroManagementDataProvider.
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
                    max_duration: {max_duration}
                TimedEvents: null
                StateEvents: null
          - {campaign_end_date}: null
        """

    def __init__(self, HDFstore, grid_no, crop_no, campaign_year, campaign_start=None):
        list.__init__(self)
        self.HDFstore = os.path.abspath(HDFstore)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)

        # Read crop name
        self.crop_name = fetch_crop_name(self.HDFstore, self.crop_no)

        # Read crop calendar
        wh = "grid_no=%i and crop_no=%i and year=%i" % (self.grid_no, self.crop_no, self.campaign_year)
        df = pd.read_hdf(self.HDFstore, "/crop_calendar", where=wh)
        if len(df) == 0:
            msg = "Failed deriving crop calendar for %s" % wh
            raise exc.PCSEError(msg)

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        start_type = str(df.start_type.iloc[0]).strip()
        self.crop_start_date = check_date(df.start_date.iloc[0])
        if start_type == "FIXED_SOWING":
            self.crop_start_type = "sowing"
        elif start_type == "FIXED_EMERGENCE":
            self.crop_start_type = "emergence"
        else:
            msg = "Unsupported START_TYPE in CROP_CALENDAR table: %s" % start_type
            raise exc.PCSEError(msg)

        # determine the campaign_start_date
        if campaign_start is None:
            self.campaign_start_date = self.crop_start_date
        elif isinstance(campaign_start, (int, float)):
            ndays = abs(int(campaign_start))
            self.campaign_start_date = self.crop_start_date - dt.timedelta(days=ndays)
        else:
            try:
                campaign_start = check_date(campaign_start)
                if campaign_start <= self.crop_start_date:
                    self.campaign_start_date = campaign_start
                else:
                    msg = "Date (%s) specified by keyword 'campaign_start' in call to AgroManagementDataProvider " \
                          "is later then crop_start_date defined in the CGMS database."
                    raise exc.PCSEError(msg % campaign_start)
            except KeyError as e:
                msg = "Value (%s) of keyword 'campaign_start' not recognized in call to AgroManagementDataProvider."
                raise exc.PCSEError(msg % campaign_start)

        # Determine crop end date/type and the end of the campaign
        self.crop_end_type = str(df.end_type.iloc[0]).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % self.crop_end_type)
            raise exc.PCSEError(msg)

        # Determine maximum duration of the crop

        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
            self.max_duration = int(df.max_duration.iloc[0])
            self.campaign_end_date = self.crop_start_date + dt.timedelta(days=self.max_duration)
        else:
            self.crop_end_date = check_date(df.end_date.iloc[0])
            self.campaign_end_date = self.crop_end_date
            self.max_duration = (self.crop_end_date - self.crop_start_date).days + 1

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
                                                     max_duration=self.max_duration,
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

        This is useful only when the INITIAL_SOIL_WATER table in CGMS12 defines a different
        campaign_start
        """
        self.campaign_start_date = check_date(start_date)
        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def __str__(self):
        msg1 = ("Agromanagement data for crop/grid/year (%s/%i/%i) derived from: %s" %
               (self.crop_name, self.grid_no, self.campaign_year, self.HDFstore))
        msg2 = self._build_yaml_agromanagement()
        msg = "  %s:\n %s" % (msg1, msg2)
        return msg


class STU_Suitability(set):
    """Returns a set() of suitable STU's for given crop_no.

    :param HDFstore: path to HDFstore to retrieve data from
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """

    def __init__(self, HDFstore, crop_no):
        self.crop_no = int(crop_no)
        self.HDFstore = os.path.abspath(HDFstore)

        # Read crop name and crop group
        wh = "crop_no=%i" % self.crop_no
        df_crop = pd.read_hdf(self.HDFstore, "/crop", where=wh)
        if len(df_crop) == 0:
            msg = "Failure reading crop type: crop_no==%i not found!" % self.crop_no
            raise exc.WeatherDataProviderError(msg)
        self.crop_name = df_crop.crop_name.iloc[0]

        self.cropgroup_no = df_crop.cropgroup_no.iloc[0]

        # read suitable STU
        wh = "cropgroup_no=%i" % self.cropgroup_no
        df = pd.read_hdf(self.HDFstore, "/suitability", where=wh)
        if len(df) == 0:
            msg = "No suitable STU found for crop_no=%s" % self.crop_no
            raise exc.PCSEError(msg)

        set.__init__(self, df.stu_no.astype(np.int))


class SiteDataProvider(dict):
    """Provides the site data from the tables INITIAL_SOIL_WATER and SITE.

    :param HDFstore: path to HDFstore to retrieve data from
    :param grid_no:  Grid number (int)
    :param crop_no: Crop number (int)
    :param campaign_year: Campaign year (int)
    :param stu_no: soil typologic unit number (int)

    Note that the parameter SSI (Initial surface storage) is
    set to zero

    Moreover, the start date of the water balance is defined by the
    column GIVEN_STARTDATE_WATBAL. This value can be accessed as
    an attribute `start_date_waterbalance`.

    """

    def __init__(self, HDFstore, grid_no, crop_no, campaign_year, stu_no):
        dict.__init__(self)

        self.HDFstore = os.path.abspath(HDFstore)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.stu_no = int(stu_no)
        self.crop_name = fetch_crop_name(self.HDFstore, self.crop_no)

        wh = "grid_no=%i and crop_no=%i and year=%i and stu_no=%i" % \
             (self.grid_no, self.crop_no, self.campaign_year, self.stu_no)
        df = pd.read_hdf(self.HDFstore, "/initial_soil_water", where=wh)
        if len(df) == 0:
            msg = "Failed retrieving site data for %s" % wh
            raise exc.PCSEError(msg)

        # Initial amount of soil water
        self["WAV"] = float(df.wav.iloc[0])

        # Raise an error in case simulation with ground water influence
        if int(df.zti.iloc[0]) != 999 or int(df.dd.iloc[0]) != 999:
            msg = ("Simulation with ground water for %s. Not implemented in PCSE/WOFOST (yet)." % wh)
            raise exc.PCSEError(msg)

        # Start date water balance
        self.start_date_waterbalance = check_date(df.given_startdate_watbal.iloc[0])

        # Derived global parameters from table SITE
        df_site = pd.read_hdf(self.HDFstore, "/site")
        self["IFUNRN"] = int(df_site.ifunrn.iloc[0])
        self["NOTINF"] = float(df_site.notinf.iloc[0])
        self["SSMAX"] = float(df_site.max_surface_storage.iloc[0])
        self["SSI"] = 0.

    def __str__(self):
        msg = ("Site parameter values for grid_no=%s, crop_no=%s (%s), stu_no=%s, "
               "campaign_year=%i derived from %s\n" % (self.grid_no, self.crop_no,
                self.crop_name, self.stu_no, self.campaign_year, self.HDFstore))
        msg += "    %s" % dict.__str__(self)

        return msg


class SoilDataIterator(list):
    """Class for iterating over the different soils in a CGMS grid.

    :param HDFstore: path to HDFstore to retrieve data from
    :param grid_no: Integer grid ID, maps to the grid_no column in the table

    Instances of this class behave like a list, allowing to iterate
    over the soils in a CGMS grid. An example::

        >>> soil_iterator = SoilDataIterator(HDFstore, grid_no=15060)
        >>> print(soildata)
        Soil data for grid_no=15060 derived from /home/wit015/test/CGMS12EU.h5
          smu_no=9050131, area=625000000, stu_no=9000282 covering 50% of smu.
            Soil parameters {'SMLIM': 0.312, 'SMFCF': 0.312, 'SMW': 0.152, 'CRAIRC': 0.06,
                             'KSUB': 10.0, 'RDMSOL': 10.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.439}
          smu_no=9050131, area=625000000, stu_no=9000283 covering 50% of smu.
            Soil parameters {'SMLIM': 0.28325, 'SMFCF': 0.28325, 'SMW': 0.12325, 'CRAIRC': 0.06,
                             'KSUB': 10.0, 'RDMSOL': 40.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.42075}
        >>> for smu_no, area, stu_no, percentage, soil_par in soildata:
        ...     print(smu_no, area, stu_no, percentage)
        ...
        (9050131, 625000000, 9000282, 50)
        (9050131, 625000000, 9000283, 50)

    """

    # name of the table with Elementary Mapping Units
    emu_table_name = "emu"

    def __init__(self, HDFstore, grid_no):

        list.__init__(self)
        self.HDFstore = os.path.abspath(HDFstore)
        self.grid_no = int(grid_no)

        SMUs = self._get_SMU_from_EMU()
        for grid_no, smu_no, area in SMUs:
            STUs = self._get_STU_from_SMU(smu_no)
            for smu_no, stu_no, percentage in STUs:
                soil_par = SoilDataProviderSingleLayer(self.HDFstore, stu_no)
                self.append((smu_no, area, stu_no, percentage, soil_par))

    def _get_SMU_from_EMU(self):
        """Retrieves the relevant SMU for given grid_no from
        table EMU.
        """
        wh = "grid_no = %i" % self.grid_no
        df_emu = pd.read_hdf(self.HDFstore, "/%s" % self.emu_table_name, where=wh)
        if len(df_emu) == 0:
            msg = ("No soil mapping units (SMU) found for grid_no=%i "
                   "in table %s" % (self.grid_no, self.emu_table_name))
            raise exc.PCSEError(msg)

        return df_emu.itertuples(index=False)

    def _get_STU_from_SMU(self, smu_no):
        """Retrieves the relevant STU for given SMU_NO from table
        SOIL_ASSOCIATION_COMPOSITION
        """
        wh = "smu_no = %i" % smu_no
        df = pd.read_hdf(self.HDFstore, "soil_association_composition", where=wh)
        if len(df) == 0:
            msg = "No soil typologic units (STU) found for smu_no=%i" % smu_no
            raise exc.PCSEError(msg)

        return df.itertuples(index=False)

    def __str__(self):
        msg = "Soil data for grid_no=%i derived from %s\n" % (self.grid_no, self.HDFstore)
        template = "  smu_no=%i, area=%.0f, stu_no=%i covering %i%% of smu.\n    Soil parameters %s\n"
        for t in self:
            msg += template % t
        return msg

