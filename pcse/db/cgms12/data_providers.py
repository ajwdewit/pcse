# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2015
"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 11 database schema.
"""

import datetime
import logging

from sqlalchemy import MetaData, select, Table, and_
from tabulate import tabulate
import numpy as np
import yaml

from ...util import wind10to2, safe_float, check_date, reference_ET
from ... import exceptions as exc
from ...base_classes import WeatherDataContainer, WeatherDataProvider


def fetch_crop_name(engine, crop_no):
    """Retrieves the name of the crop from the CROP table for
    given crop_no.

    :param engine: SqlAlchemy engine object providing DB access
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """
    metadata = MetaData(engine)
    table_crop = Table("crop", metadata, autoload=True)
    r = select([table_crop],
               table_crop.c.crop_no == crop_no).execute()
    row = r.fetchone()
    r.close()
    if row is None:
        msg = "Failed deriving crop name from CROP table for crop_no %s" % crop_no
        raise exc.PCSEError(msg)
    return row.crop_name


class STU_Suitability(set):
    """Returns a set() of suitable STU's for given crop_no.

    :param engine: SqlAlchemy engine object providing DB access
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """

    def __init__(self, engine, crop_no):
        self.crop_no = int(crop_no)

        metadata = MetaData(engine)
        table_crop = Table('crop', metadata, autoload=True)
        r = select([table_crop], table_crop.c.crop_no == self.crop_no).execute()
        row = r.fetchone()
        if row is None:
            msg = "No crop group found for crop_no=%s" % self.crop_no
            raise exc.PCSEError(msg)
        self.cropgroup_no = int(row.cropgroup_no)
        self.crop_name = row.crop_name

        table_ss = Table('suitability', metadata, autoload=True, )
        r = select([table_ss.c.stu_no], table_ss.c.cropgroup_no == self.cropgroup_no).execute()
        rows = r.fetchall()
        set.__init__(self, [int(row.stu_no) for row in rows])


class WeatherObsGridDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the WEATHER_OBS_GRID table in a CGMS12
    compatible database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no:  Grid number (int) to retrieve data for
    :param start_date: Retrieve meteo data starting with start_date
        (datetime.date object)
    :param end_date: Retrieve meteo data up to and including end_date
        (datetime.date object)
    :param recalc_ET: Set to True to force calculation of reference
        ET values. Mostly useful when values have not been calculated
        in the CGMS database.

    Note that all meteodata is first retrieved from the DB and stored
    internally. Therefore, no DB connections are stored within the class
    instance. This makes that class instances can be pickled.

    If start_date and end_date are not provided then the entire time-series
    for the grid is retrieved.
    """
    # default values for the Angstrom parameters in the sunshine duration model
    angstA = 0.18
    angstB = 0.55
    def __init__(self, engine, grid_no, start_date=None, end_date=None,
                 recalc_ET=False):

        WeatherDataProvider.__init__(self)

        self.grid_no = grid_no
        self.recalc_ET = recalc_ET

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
        """Retrieves latitude, longitude, elevation from "grid" table and
        assigns them to self.latitude, self.longitude, self.elevation."""

        # Pull Latitude value for grid nr from database

        table_grid = Table("grid", metadata, autoload=True)
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
        """Retrieves the meteo data from table "grid_weather".
        """

        try:
            # if start_date/end_date are None, define a date in the far past/future
            start_date = self.start_date if self.start_date is not None else datetime.date(1, 1, 1)
            end_date = self.end_date if self.end_date is not None else datetime.date(9999, 1, 1)
            table_db = Table("weather_obs_grid", metadata, autoload=True)
            r = select([table_db], and_(table_db.c.grid_no == self.grid_no,
                                        table_db.c.day >= start_date,
                                        table_db.c.day <= end_date)
                       ).execute()
            rows = r.fetchall()

            c = len(rows)
            if self.time_interval is not None:
                if c < self.time_interval:
                    msg = ("Only %i records selected from table 'WEATHER_OBS_GRID' "
                           "for grid %i, period %s -- %s.")
                    self.logger.warn(msg, c, self.grid_no, self.start_date,
                                     self.end_date)

            for row in rows:
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = self._make_WeatherDataContainer(row, t)
                self._store_WeatherDataContainer(wdc, DAY)
        except Exception:
            msg = "Failure reading meteodata for grid %s "
            self.logger.exception(msg, self.grid_no)
            raise exc.PCSEError(msg, self.grid_no)

        msg = ("Successfully retrieved weather data from 'WEATHER_OBS_GRID' table "
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
                  "IRRAD": float(row.radiation)*1000.,
                  "SNOWDEPTH": safe_float(row.snowdepth)})

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
    """Class for providing agromanagement data from the CROP_CALENDAR table in a CGMS12 database.

    :param engine: SqlAlchemy engine object providing DB access
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
                    crop_id: '{crop_name}'
                    crop_start_date: {crop_start_date}
                    crop_start_type: {crop_start_type}
                    crop_end_date: {crop_end_date}
                    crop_end_type: {crop_end_type}
                    max_duration: {max_duration}
                TimedEvents: null
                StateEvents: null
          - {campaign_end_date}: null
        """

    def __init__(self, engine, grid_no, crop_no, campaign_year, campaign_start=None):
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

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        cgms12_start_type = str(row.start_type).strip()
        self.crop_start_date = check_date(row.start_date)
        if cgms12_start_type == "FIXED_SOWING":
            self.crop_start_type = "sowing"
        elif cgms12_start_type == "FIXED_EMERGENCE":
            self.crop_start_type = "emergence"
        else:
            msg = "Unsupported START_TYPE in CROP_CALENDAR table: %s" % row.start_type
            raise exc.PCSEError(msg)

        # determine the campaign_start_date
        if campaign_start is None:
            self.campaign_start_date = self.crop_start_date
        elif isinstance(campaign_start, (int, float)):
            ndays = abs(int(campaign_start))
            self.campaign_start_date = self.crop_start_date - datetime.timedelta(days=ndays)
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
        self.crop_end_type = str(row.end_type).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % row.end_type)
            raise exc.PCSEError(msg)

        # Determine maximum duration of the crop

        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
            self.max_duration = int(row.max_duration)
            self.campaign_end_date = self.crop_start_date + datetime.timedelta(days=self.max_duration)
        else:
            self.crop_end_date = check_date(row.end_date)
            self.campaign_end_date = self.crop_end_date
            self.max_duration = (self.crop_end_date - self.crop_start_date).days + 1

        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def _build_yaml_agromanagement(self):
        """Builds the YAML agromanagent string"""

        input = self.agro_management_template.format(campaign_start_date=self.campaign_start_date,
                                                     crop_name=self.crop_name,
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
            items = yaml.load(input)
        except yaml.YAMLError as e:
            msg = "Failed parsing agromanagement string %s: %s" % (input, e)
            raise exc.PCSEError(msg)
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
        msg1 = ("Agromanagement data for crop_no=%i (%s) derived from: %s" %
               (self.crop_no, self.crop_name, self.db_resource))
        msg2 = self._build_yaml_agromanagement()
        msg = "  %s:\n %s" % (msg1, msg2)
        return msg


class SoilDataProviderSingleLayer(dict):
    """Class for providing soil data from the ROOTING_DEPTH AND
    SOIL_PHYSICAL_GROUP tableS in a CGMS8/12 database. This
    applies to the single layered soil only.

    :param engine: SqlAlchemy engine object providing DB access
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

    def __init__(self, engine, stu_no):
        dict.__init__(self)

        metadata = MetaData(engine)

        # First get the rooting depth class and soil_group_no
        rd_class, spg_no = self._get_from_STU(metadata, stu_no)
        # Get the actual rooting depth [cm]
        self._get_rooting_depth(metadata, rd_class)
        # Get the actual soil hydrological parameters.
        self._get_soil_hydraulic_parameters(metadata, spg_no)

        # define SMLIM
        self["SMLIM"] = self["SMFCF"]

    def _get_from_STU(self, metadata, stu_no):
        """Retrieves the soil parameters for the given soil typologic unit
        (stu_no) from the tables SOIL_PHYSICAL_GROUP and ROOTING_DEPTH.
        """
        table_stu = Table("soil_typologic_unit", metadata, autoload=True)

        # first get the soil_group_no and rooting depth class
        s = select([table_stu.c.soil_group_no, table_stu.c.calculated_rooting_depth],
                   table_stu.c.stu_no == stu_no).execute()
        row = s.fetchone()
        if row is None:
            msg = ("No record found for stu_no=%i in table "
                   "SOIL_TYPOLOGIC_UNIT." % stu_no)
            raise exc.PCSEError(msg)
        soil_group_no = int(row.soil_group_no)
        rd_class = int(row.calculated_rooting_depth)

        return (rd_class, soil_group_no)

    def _get_rooting_depth(self, metadata, rd_class):
        """Gets the rooting depth from the table ROOTING_DEPTH and
         stores into self[] directly under parameter name 'RDMSOL'.

        :param metadata: An SQLAlchemy Metadata object
        :param rd_class: The rooting depth class (integer)
        """

        table_rd = Table("rooting_depth", metadata, autoload=True)
        s = select([table_rd]).execute()
        rows = s.fetchall()
        # Note that we need to loop over the row here instead of putting a
        # WHERE class = rd_class in the SQL query because 'class' is a reserved
        # word and the column expression table_rd.c.class raises SyntaxError
        for c, root_depth in rows:
            if rd_class == c:
                break
        else:
            msg = ("No record found for rooting_depth_class=%i in table "
                   "ROOTING_DEPTH." % rd_class)
            raise exc.PCSEError(msg)

        self["RDMSOL"] = float(root_depth)

    def _get_soil_hydraulic_parameters(self, metadata, spg_no):
        """Retrieves the soil hydraulic parameters and stores into self[] directly.

        :param metadata: An SQLAlchemy Metadata object
        :param spg_no: the soil physical group number (integer)
        :return: None
        """
        table_spg = Table("soil_physical_group", metadata, autoload=True)

        for (wofost_soil_par, db_soil_par) in self.soil_parameters:
            r = select([table_spg],
                       and_(table_spg.c.soil_group_no == spg_no,
                            table_spg.c.parameter_code == db_soil_par)).execute()
            row = r.fetchone()
            if row is None:
                msg = "Parameter %s not found in table SOIL_PHYSICAL_GROUP" % db_soil_par
                raise exc.PCSEError(msg)
            self[wofost_soil_par] = float(row.parameter_xvalue)


class SoilDataIterator(list):
    """Class for iterating over the different soils in a CGMS grid.

    Instances of this class behave like a list, allowing to iterate
    over the soils in a CGMS grid. An example::

        >>> soil_iterator = SoilDataIterator(engine, grid_no=15060)
        >>> print(soildata)
        Soil data for grid_no=15060 derived from oracle+cx_oracle://cgms12eu:***@eurdas.world
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

    def __init__(self, engine, grid_no):

        list.__init__(self)
        metadata = MetaData(engine)
        self.db_resource = str(engine)[7:-1]
        self.grid_no = int(grid_no)

        SMUs = self._get_SMU_from_EMU(metadata, self.grid_no)
        for smu_no, area in SMUs:
            STUs = self._get_STU_from_SMU(metadata, smu_no)
            for stu_no, percentage in STUs:
                soil_par = SoilDataProviderSingleLayer(engine, stu_no)
                self.append((smu_no, area, stu_no, percentage, soil_par))

    def _get_SMU_from_EMU(self, metadata, grid_no):
        """Retrieves the relevant SMU for given grid_no from
        table EMU.
        """
        table_emu = Table(self.emu_table_name, metadata, autoload=True)
        r = select([table_emu.c.smu_no, table_emu.c.area],
                   table_emu.c.grid_no == grid_no).execute()
        rows = r.fetchall()
        if rows is None:
            msg = ("No soil mapping units (SMU) found for grid_no=%i "
                   "in table %s" % (grid_no, self.emu_table_name))
            raise exc.PCSEError(msg)

        return rows

    def _get_STU_from_SMU(self, metadata, smu_no):
        """Retrieves the relevant STU for given SMU_NO from table
        SOIL_ASSOCIATION_COMPOSITION
        """
        table_sac = Table("soil_association_composition", metadata, autoload=True)
        r = select([table_sac.c.stu_no, table_sac.c.percentage],
                   table_sac.c.smu_no == smu_no).execute()
        rows = r.fetchall()
        if rows is None:
            msg = "No soil typologic units (STU) found for smu_no=%i" % smu_no
            raise exc.PCSEError(msg)

        return rows

    def __str__(self):
        msg = "Soil data for grid_no=%i derived from %s\n" % (self.grid_no, self.db_resource)
        template = "  smu_no=%i, area=%.0f, stu_no=%i covering %i%% of smu.\n    Soil parameters %s\n"
        for t in self:
            msg += template % t
        return msg

class CropDataProvider(dict):
    """Retrieves the crop parameters for the given grid_no, crop_no and year
    from the tables CROP_CALENDAR, CROP_PARAMETER_VALUE and VARIETY_PARAMETER_VALUE.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the GRID_NO column in the table
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """

    # Define single and tabular crop parameter values
    parameter_codes_single = ("CFET", "CVL", "CVO", "CVR", "CVS", "DEPNR", "DLC",
                              "DLO", "DVSEND", "EFF", "IAIRDU", "IDSL", "KDIF",
                              "LAIEM", "PERDL", "Q10", "RDI", "RDMCR", "RGRLAI",
                              "RML", "RMO", "RMR", "RMS", "RRI", "SPA", "SPAN", "SSA",
                              "TBASE", "TBASEM", "TDWI", "TEFFMX", "TSUM1", "TSUM2",
                              "TSUMEM")
    parameter_codes_tabular = ("AMAXTB", "DTSMTB", "FLTB", "FOTB", "FRTB", "FSTB",
                               "RDRRTB", "RDRSTB", "RFSETB", "SLATB", "TMNFTB",
                               "TMPFTB")
    # Some parameters have to be converted from a single to a tabular form
    single2tabular = {"SSA": ("SSATB", [0., None, 2.0, None]),
                      "KDIF": ("KDIFTB", [0., None, 2.0, None]),
                      "EFF": ("EFFTB", [0., None, 40., None])}
    # Default values for additional parameters not defined in CGMS
    parameters_additional = {"DVSI": 0.0, "IOX": 0}

    def __init__(self, engine, grid_no, crop_no, campaign_year):
        dict.__init__(self)

        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.crop_name = fetch_crop_name(engine, self.crop_no)
        self.db_resource = str(engine)[7:-1]

        metadata = MetaData(engine)

        # Get crop variety from crop_calendar;
        table_cc = Table('crop_calendar', metadata, autoload=True)
        r = select([table_cc],
                   and_(table_cc.c.grid_no == self.grid_no,
                        table_cc.c.crop_no == self.crop_no,
                        table_cc.c.year == self.campaign_year)).execute()
        row = r.fetchone()
        if row is None:
            msg = "No entry found in CROP_CALENDAR for grid_no=%s, crop_no=%s, year=%s."
            raise exc.PCSEError(msg % (self.grid_no, self.crop_no, self.campaign_year))
        self.variety_no = int(row.variety_no)

        # get parameters from CGMS db
        self._fetch_crop_parameter_values(metadata, self.crop_no)
        self._fetch_variety_parameter_values(metadata, self.crop_no, self.variety_no)
        self.update(self.parameters_additional)

        # Finally add crop name
        self["CRPNAM"] = self.crop_name

    def _fetch_crop_parameter_values(self, metadata, crop_no):
        """Derived the crop parameter values from the CROP_PARAMETER_VALUE
        table for given crop_no and add directly to dict self[]..
        """

        # Pull single value parameters from CROP_PARAMETER_VALUE
        table_crop_pv = Table('crop_parameter_value', metadata, autoload=True)
        for parameter_code in self.parameter_codes_single:
            r = select([table_crop_pv],
                       and_(table_crop_pv.c.crop_no == crop_no,
                            table_crop_pv.c.parameter_code == parameter_code)).execute()
            row = r.fetchone()
            if row is None:
                msg = "No parameter value found for crop_no=%s, parameter_code='%s'."
                raise exc.PCSEError(msg % (self.crop_no, parameter_code))
            if parameter_code not in self.single2tabular:
                self[parameter_code] = float(row.parameter_xvalue)
            else:
                pvalue = float(row.parameter_xvalue)
                code, value = self._convert_single2tabular(parameter_code, pvalue)
                self[code] = value

        # Pull tabular parameters from CROP_PARAMETER_VALUE
        for parameter_code in self.parameter_codes_tabular:
            pattern = parameter_code + r'%'
            r = select([table_crop_pv],
                       and_(table_crop_pv.c.crop_no == crop_no,
                            table_crop_pv.c.parameter_code.like(pattern)),
                       order_by=[table_crop_pv.c.parameter_code]).execute()
            rows = r.fetchall()
            if not rows:
                msg = "No parameter value found for crop_no=%s, parameter_code='%s'."
                raise exc.PCSEError(msg % (self.crop_no, parameter_code))
            if len(rows) == 1:
                msg = ("Single parameter value found for crop_no=%s, parameter_code='%s' while "
                       "tabular parameter expected." % (crop_no, parameter_code))
                raise exc.PCSEError(msg)
            values = []
            for row in rows:
                values.extend([float(row.parameter_xvalue), float(row.parameter_yvalue)])
            self[parameter_code] = values

    def _fetch_variety_parameter_values(self, metadata, crop_no, variety_no):
        """Derived the crop parameter values from the VARIETY_PARAMETER_VALUE
        table for given crop_no & variety_on and add directly to dict self[].
        """

        # Pull single value parameters from VARIETY_PARAMETER_VALUE
        table_crop_vpv = Table('variety_parameter_value', metadata, autoload=True)
        for parameter_code in self.parameter_codes_single:
            r = select([table_crop_vpv],
                       and_(table_crop_vpv.c.crop_no == crop_no,
                            table_crop_vpv.c.variety_no == variety_no,
                            table_crop_vpv.c.parameter_code == parameter_code)).execute()
            row = r.fetchone()
            if row is None:
                continue

            if parameter_code not in self.single2tabular:
                self[parameter_code] = float(row.parameter_xvalue)
            else:
                pvalue = float(row.parameter_xvalue)
                code, value = self._convert_single2tabular(parameter_code, pvalue)
                self[code] = value

        # Pull tabular parameters from CROP_PARAMETER_VALUE
        for parameter_code in self.parameter_codes_tabular:
            pattern = parameter_code + r'%'
            r = select([table_crop_vpv],
                       and_(table_crop_vpv.c.crop_no == crop_no,
                            table_crop_vpv.c.variety_no == variety_no,
                            table_crop_vpv.c.parameter_code.like(pattern)),
                       order_by=[table_crop_vpv.c.parameter_code]).execute()
            rows = r.fetchall()
            if not rows:
                continue

            if len(rows) == 1:
                msg = ("Single parameter value found for crop_no=%s, parameter_code='%s' while "
                       "tabular parameter expected." % (crop_no, parameter_code))
                raise exc.PCSEError(msg)
            values = []
            for row in rows:
                values.extend([float(row.parameter_xvalue), float(row.parameter_yvalue)])
            self[parameter_code] = values

    def _convert_single2tabular(self, parameter_code, pvalue):
        """Converts the single parameter into a tabular parameter.
        """
        tabular_parameter_code, template = self.single2tabular[parameter_code]
        tabular_values = [pvalue if v is None else v for v in template]

        return tabular_parameter_code, tabular_values

    def __str__(self):
        msg = ("Crop parameter values for grid_no=%s, crop_no=%s (%s), variety_no=%s, "
               "campaign_year=%i derived from %s\n" %
               (self.grid_no, self.crop_no, self.crop_name, self.variety_no,
                self.campaign_year, self.db_resource))
        single_values = []
        tabular_values = []
        for pcode in sorted(self.keys()):
            value = self[pcode]
            if not isinstance(value, list):
                single_values.append((pcode, value))
            else:
                tabular_values.append((pcode, value))

        # Format the single parameters in a table of 4 columns
        msg += "Single parameter values:\n"
        # If not of even length add ["",""]
        if not len(single_values) % 2 == 0:
            single_values.append(["", ""])
        np_single_values = np.array(single_values, dtype=np.string_)
        shp = np_single_values.shape
        np_single_values.shape = (shp[0]/2, shp[1]*2)
        msg += tabulate(np_single_values, headers=["Par_code", "Value", "Par_code", "Value"])
        msg += "\n"

        # Format the tabular parameters in two columns
        msg += "Tabular parameters:\n"
        msg += tabulate(tabular_values, headers=["Par_code", "Value"])
        msg += "\n"

        return msg


class SiteDataProvider(dict):
    """Provides the site data from the tables INITIAL_SOIL_WATER and SITE.

    :param engine: SqlAlchemy engine object providing DB access
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

    def __init__(self, engine, grid_no, crop_no, campaign_year, stu_no):
        dict.__init__(self)

        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.stu_no = int(stu_no)
        self.crop_name = fetch_crop_name(engine, self.crop_no)
        self.db_resource = str(engine)[7:-1]

        metadata = MetaData(engine)
        table_isw = Table('initial_soil_water', metadata, autoload=True)
        r = select([table_isw], and_(table_isw.c.grid_no == self.grid_no,
                                     table_isw.c.crop_no == self.crop_no,
                                     table_isw.c.year == self.campaign_year,
                                     table_isw.c.stu_no == self.stu_no)).execute()
        row = r.fetchone()
        r.close()
        if row is None:
            msg = ("Failed retrieving site data for grid_no=%s, crop_no=%s, "
                   "campaign_year=%s, stu_no=%s" % (self.grid_no, self.crop_no,
                                                    self.campaign_year, self.stu_no))
            raise exc.PCSEError(msg)
        # Initial amount of soil water
        self["WAV"] = float(row.wav)

        # Raise an error in case simulation with ground water influence
        if int(row.zti) != 999 or int(row.dd) != 999:
            msg = ("Simulation with ground water for grid_no=%s, crop_no=%s, "
                   "campaign_year=%s, stu_no=%s. Not implemented in PCSE/WOFOST (yet)."
                   % (self.grid_no, self.crop_no, self.campaign_year, self.stu_no))
            raise exc.PCSEError(msg)

        # Start date water balance
        self.start_date_waterbalance = check_date(row.given_startdate_watbal)

        # Derived global parameters from table SITE
        table_site = Table('site', metadata, autoload=True)
        r = select([table_site]).execute()
        row = r.fetchone()
        self["IFUNRN"] = int(row.ifunrn)
        self["NOTINF"] = int(row.notinf)
        self["SSMAX"] = float(row.max_surface_storage)
        self["SSI"] = 0.

    def __str__(self):
        msg = ("Site parameter values for grid_no=%s, crop_no=%s (%s), stu_no=%s, "
               "campaign_year=%i derived from %s\n" % (self.grid_no, self.crop_no,
                self.crop_name, self.stu_no, self.campaign_year, self.db_resource))
        msg += "    %s" % dict.__str__(self)

        return msg
