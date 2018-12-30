# -*- coding: utf-8 -*-
# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2016

"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 14 database schema.
"""
import datetime as dt
import os

import numpy as np
import yaml
from sqlalchemy import MetaData, select, Table, and_

from ... import exceptions as exc
from ...base import WeatherDataContainer, WeatherDataProvider
from ...util import wind10to2, safe_float, check_date, reference_ET
from ... import settings
from .. import wofost_parameters


def fetch_crop_name(engine, idcrop_parametrization):
    """Retrieves the name of the crop from the crop_parametrizations table for a given idcrop_parametrization.

    :param engine: SqlAlchemy engine object providing DB access
    :param idcrop_parametrization: Integer id of crop parametrization, maps to the like-named column in the table
    """

    metadata = MetaData(engine)
    table_crop = Table("crop_parametrizations", metadata, autoload=True)
    sm = select([table_crop], table_crop.c.idcrop_parametrization == idcrop_parametrization)
    sc = sm.execute()
    row = sc.fetchone()
    sc.close()
    if row is None:
        msg = "Failed deriving crop name from view link_crop_parametrization for " \
              "idcrop_parametrization %s" % idcrop_parametrization
        raise exc.PCSEError(msg)
    result = row.crop_parametrization
    return result


class STU_Suitability(set):
    """Returns a set() of suitable STU's for a given idcrop_parametrization.

    :param engine: SqlAlchemy engine object providing DB access
    :param crop_parametrization: Integer id of crop, maps to the like-named column in the table
    """

    def __init__(self, engine, idcrop_parametrization):
        # Initialise
        self.idcrop_parametrization = int(idcrop_parametrization)
        self.crop_name = fetch_crop_name(engine, idcrop_parametrization)
        metadata = MetaData(engine)

        # Retrieve suitable STU's - first connect to the 2 relevant tables      
        t1 = Table("link_crop_aggregation", metadata, autoload=True)
        t2 = Table("link_crop_stu", metadata, autoload=True)

        # Assume that records in table link_crop_aggregations with higher idcrop are less
        # relevant than those with higher values
        sm = select([t1.c.idcrop, t1.c.idcrop_parametrization, t2],
                    and_(t1.c.idcrop_parametrization == idcrop_parametrization,
                         t1.c.idcrop == t2.c.idcrop)).order_by(t1.c.idcrop)
        sc = sm.execute()
        rows = sc.fetchall()
        sc.close()
        if not rows:
            msg = "No suitable soil type unit found for idcrop_parametrization=%s"
            raise exc.PCSEError(msg % idcrop_parametrization)
        set.__init__(self, [int(row.idstu) for row in rows])


class WeatherObsGridDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the WEATHER_ERA_GRID table in a CGMS14
    compatible database - or from another table when specified

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid:  Grid number (int) to retrieve data for
    :param start_date: Retrieve meteo data starting with start_date
        (datetime.date object)
    :param end_date: Retrieve meteo data up to and including end_date
        (datetime.date object)
    :param recalc_ET: Set to True to force calculation of reference
        ET values. Mostly useful when values have not been calculated
         in the CGMS database.
    :param use_cache: Set to False to ignore reading and writing a cache file
    :param table_name:
    """
    # default values for the Angstrom parameters in the sunshine duration model
    angstA = 0.29
    angstB = 0.49

    def __init__(self, engine, idgrid, start_date=None, end_date=None, use_cache=True,
                 recalc_ET=False, recalc_TEMP=False, table_name='weather_era_grid'):
        # Initialise
        WeatherDataProvider.__init__(self)
        self.idgrid = idgrid
        self.recalc_ET = recalc_ET
        self.recalc_TEMP = recalc_TEMP
        self.table_name = table_name
        self.use_cache = use_cache

        if not self._self_load_cache(self.idgrid) or self.use_cache is False:
            metadata = MetaData(engine)

            # Check the start and end dates and assign when ok
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

            # Get location info (lat/lon/elevation)
            self._fetch_location_from_db(metadata)

            # Retrieve the meteo data
            self._fetch_weather_from_db(metadata)

            # Provide a description that is shown when doing a print()
            line1 = "Weather data retrieved from CGMS 14 db %s" % str(engine)[7:-1]
            line2 = "for idgrid: %s" % self.idgrid
            self.description = [line1, line2]

            # Save cache file
            if self.use_cache:
                fname = self._get_cache_filename(self.idgrid)
                self._dump(fname)

    def _get_cache_filename(self, grid_no):
        fname = "%s_grid_%i.cache" % (self.__class__.__name__, grid_no)
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, fname)
        return cache_filename

    def _self_load_cache(self, grid_no):
        """Checks if a cache file exists and tries to load it."""
        cache_fname = self._get_cache_filename(grid_no)
        if os.path.exists(cache_fname):
            r = os.stat(cache_fname)
            cache_file_date = dt.date.fromtimestamp(r.st_mtime)
            age = (dt.date.today() - cache_file_date).days
            if age < 1:
                try:
                    self._load(cache_fname)
                    return True
                except exc.PCSEError:
                    pass
        return False

    # ---------------------------------------------------------------------------
    def _fetch_location_from_db(self, metadata):
        """Retrieves latitude, longitude, elevation from "grid" table and
        assigns them to self.latitude, self.longitude, self.elevation."""
        tg = Table("grids", metadata, autoload=True)
        sc = select([tg.c.latitude, tg.c.longitude, tg.c.altitude],
                    tg.c.idgrid == self.idgrid).execute()
        row = sc.fetchone()
        sc.close()
        if row is None:
            msg = "Failed deriving location info for grid %s" % self.idgrid
            raise exc.PCSEError(msg)

        # Use the resulta
        self.latitude = float(row.latitude)
        self.longitude = float(row.longitude)
        self.elevation = float(row.altitude)

        # Report success
        msg = ("Succesfully retrieved location information from 'grids' table for grid %s")
        self.logger.info(msg, self.idgrid)

    # ---------------------------------------------------------------------------
    def _fetch_weather_from_db(self, metadata):
        """Retrieves the meteo data from table "grid_weather".
        """
        try:
            # Check input - if start_date/end_date are None, define a date in the far past/future
            start_date = self.start_date if self.start_date is not None else dt.date(dt.MINYEAR, 1, 1)
            end_date = self.end_date if self.end_date is not None else dt.date(dt.MAXYEAR, 1, 1)

            # Get hold of weather table, statement, sqlalchemy cursor and finally of the rows
            tw = Table(self.table_name, metadata, autoload=True)
            sm = select([tw], and_(tw.c.idgrid == self.idgrid,
                                   tw.c.day >= start_date,
                                   tw.c.day <= end_date))
            sc = sm.execute()
            rows = sc.fetchall()

            count = len(rows)
            if self.time_interval is not None:
                if count < self.time_interval:
                    msg = ("Only %i records selected from table 'WEATHER_OBS_GRID' "
                           "for grid %i, period %s -- %s.")
                    self.logger.warn(msg, count, self.idgrid, self.start_date,
                                     self.end_date)

            for row in rows:
                DAY = self.check_keydate(row.day)
                t = {"DAY": DAY, "LAT": self.latitude,
                     "LON": self.longitude, "ELEV": self.elevation}
                wdc = self._make_WeatherDataContainer(row, t)
                self._store_WeatherDataContainer(wdc, DAY)

        except Exception:
            msg = "Failure reading meteodata for grid %s "
            self.logger.exception(msg, self.idgrid)
            raise exc.PCSEError(msg, self.idgrid)

        # Report success
        msg = "Successfully retrieved weather data from table '%s' for grid %s between %s and %s"
        self.logger.info(msg, self.table_name, self.idgrid, self.start_date, self.end_date)

    # ---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""
        result = None

        t.update({"TMAX": float(row.temperature_max),
                  "TMIN": float(row.temperature_min),
                  "TEMP": float(row.temperature_avg),
                  "VAP": float(row.vapourpressure),
                  "WIND": wind10to2(float(row.windspeed)),
                  "RAIN": float(row.precipitation) / 10.,
                  "IRRAD": float(row.radiation) * 1000.,
                  "SNOWDEPTH": safe_float(row.snowdepth)})

        if not self.recalc_ET:
            t.update({"E0": float(row.e0) / 10.,
                      "ES0": float(row.es0) / 10.,
                      "ET0": float(row.et0) / 10.})
        else:
            e0, es0, et0 = reference_ET(ANGSTA=self.angstA,
                                        ANGSTB=self.angstB, **t)

            t.update({"E0": e0 / 10.,
                      "ES0": es0 / 10.,
                      "ET0": et0 / 10.})

        if self.recalc_TEMP:
            t["TEMP"] = (float(row.temperature_max) + float(row.temperature_min))/2.

        result = WeatherDataContainer(**t)
        return result


class AgroManagementDataProvider(list):
    """Class for providing agromanagement from the CROP_CALENDAR table in a CGMS14 database.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid: Integer grid ID, maps to the idgrid column in the table
    :param idcrop_prmtrz: Integer id of crop parametrization, maps to the idcrop_parameterization
           column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """
    agro_management_template = """
          - {0[campaign_start_date]}:
                CropCalendar:
                    crop_name: '{0[crop_name]}'
                    variety_name: '{0[variety_name]}'
                    crop_start_date: {0[crop_start_date]}
                    crop_start_type: {0[start_period]}
                    crop_end_date: {0[crop_end_date]}
                    crop_end_type: {0[end_period]}
                    max_duration: {0[duration]}
                TimedEvents: null
                StateEvents: null
        """

    def __init__(self, engine, idgrid, idcrop_parametrization, campaign_year, campaign_start=None):
        # Initialise
        list.__init__(self)
        self.idgrid = idgrid
        self.idcrop_parametrization = idcrop_parametrization
        self.crop_name = fetch_crop_name(engine, idcrop_parametrization)
        self.campaign_year = campaign_year
        self.amdict = {}

        # Use the idcrop_parametrization to search in the table crop_calendars 
        metadata = MetaData(engine)
        t = Table("crop_calendars", metadata, autoload=True)
        sm = select([t], and_(t.c.idgrid == self.idgrid,
                              t.c.idcrop_parametrization == self.idcrop_parametrization,
                              t.c.year == self.campaign_year))
        sc = sm.execute()
        row = sc.fetchone()
        sc.close()

        # Process the query result - dates should be in the format 'yyyy-mm-dd'!
        if row is None:
            msg = "Failed deriving agromanagement info for grid %s" % self.idgrid
            raise exc.PCSEError(msg)

        for key, value in row.items():
            if value in ["EMERGENCE", "SOWING", "HARVEST", "MATURITY"]:
                value = value.lower()
            if key == "duration":
                value = int(value)
            self.amdict[key] = value

        self.conditional_datecopy("start_period", "crop_start_date", "emergence", "sowing")
        self.conditional_datecopy("end_period", "crop_end_date", "maturity", "harvesting")
        self.amdict["campaign_start_date"] = check_date(self.amdict["crop_start_date"])
        self.amdict["campaign_end_date"] = check_date(self.amdict["crop_end_date"]) + dt.timedelta(days=1)
        self.amdict["crop_name"] = self.crop_name
        # We do not get a variety_name from the CGMS database, so we make one
        # as <crop_name>_<grid>_<year>
        self.amdict["variety_name"] = "%s_%s_%s" % (self.crop_name, self.idgrid, self.campaign_year)

        # determine the campaign_start_date
        if campaign_start is None:
            self.amdict["campaign_start_date"] = self.amdict['crop_start_date']
        elif isinstance(campaign_start, (int, float)):
            ndays = abs(int(campaign_start))
            self.amdict["campaign_start_date"] = self.amdict["crop_start_date"] - dt.timedelta(days=ndays)
        else:
            try:
                campaign_start = check_date(campaign_start)
                if campaign_start <= self.amdict["crop_start_date"]:
                    self.amdict["campaign_start_date"] = campaign_start
                else:
                    msg = "Date (%s) specified by keyword 'campaign_start' in call to AgroManagementDataProvider " \
                          "is later then crop_start_date defined in the CGMS database."
                    raise exc.PCSEError(msg % campaign_start)
            except KeyError as e:
                msg = "Value (%s) of keyword 'campaign_start' not recognized in call to AgroManagementDataProvider."
                raise exc.PCSEError(msg % campaign_start)

        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def conditional_datecopy(self, testkey, targetkey, value1, value2):
        if self.amdict == None: return
        try:
            if self.amdict[testkey].lower() == value1.lower():
                self.amdict[targetkey] = check_date(self.amdict[value1])
            else:
                self.amdict[targetkey] = check_date(self.amdict[value2])
        except KeyError:
            raise exc.PCSEError("Failed to add value to dictionary for key %s" % targetkey)

    def _build_yaml_agromanagement(self):
        """Builds the YAML agromanagent string"""

        return self.agro_management_template.format(self.amdict)

    def _parse_yaml(self, input):
        """Parses the input YAML string and assigns to self"""
        try:
            items = yaml.load(input)
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
        self.amdict["campaign_start_date"] = check_date(start_date)
        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)


class SoilDataProviderSingleLayer(dict):
    """Class for providing soil data from the ROOING_DEPTH AND
    SOIL_PHYSICAL_GROUP tableS in a CGMS14 database. This
    applies to the single layered soil only.

    :param engine: SqlAlchemy engine object providing DB access
    :param idstu: Integer idstu, maps to the idstu column in the table
                   link_weighted_parameters for providing soil data

    Note that the value of the parameter SMLIM (Initial maximum
    moisture content in initial rooting depth zone) is set
    to field capacity (SMFCF)
    """
    soil_hydraulic_parameters = [("CRAIRC", "CRITICAL_AIR_CONTENT"),
                                 ("K0", "HYDR_CONDUCT_SATUR"),
                                 ("SOPE", "MAX_PERCOL_ROOT_ZONE"),
                                 ("KSUB", "MAX_PERCOL_SUBSOIL")]
    soil_moisture_content_parameters = [("SMFCF", "soil_moisture_fc"),
                                        ("SM0", "soil_moisture_sat"),
                                        ("SMW", "soil_moisture_wp"),
                                        ("RDMSOL", "depth")]  # CALCULATED_ROOTING_?

    def __init__(self, engine, idstu):
        # Initialise
        dict.__init__(self)
        metadata = MetaData(engine)

        # Get the actual rooting depth [cm]
        self._soil_moisture_content_parameters(metadata, idstu)
        # Get the actual soil hydrological parameters.
        self._get_soil_hydraulic_parameters(metadata, idstu)

        # define SMLIM
        self["SMLIM"] = self["SMFCF"]

    def _soil_moisture_content_parameters(self, metadata, idstu):
        """Gets the soil moisture content parameters from the table 
         link_weighted_parameters and them stores into self[] directly.

        :param metadata: An SQLAlchemy Metadata object
        :param idstu: the id for the soil typologic unit(integer)
        """
        t = Table("link_stu_weighted_parameters", metadata, autoload=True)
        sc = select([t], t.c.idstu == idstu).execute()
        row = sc.fetchone()
        sc.close()
        if row is None:
            msg = "No soil moisture content parameters found in table " \
                  "link_stu_weighted_parameters for idstu=%s" % idstu
            raise exc.PCSEError(msg)

        for (wofost_soil_par, db_soil_par) in self.soil_moisture_content_parameters:
            self[wofost_soil_par] = float(getattr(row, db_soil_par))

    def _get_soil_hydraulic_parameters(self, metadata, idstu):
        """Defineds the soil hydraulic parameters and stores into self[] directly.

        :param metadata: An SQLAlchemy Metadata object
        :param idstu: the soil physical group number (integer)
        :return: None

        NOTE: soil hydraulic parameters are not defined in BioMA database yet.
              therefore, default values are used.
        """
        for (wofost_soil_par, db_soil_par) in self.soil_hydraulic_parameters:
            if wofost_soil_par == "CRAIRC":
                self[wofost_soil_par] = 0.06
            else:
                self[wofost_soil_par] = 10.0


class SoilDataIterator(list):
    """Class for iterating over the different soils in a CGMS grid.

    Instances of this class behave like a list, allowing to iterate
    over the soils in a CGMS grid. An example::

    >>> soil_iterator = SoilDataIterator(engine, idgrid=15060)
    >>> print(soildata)
    Soil data for idgrid=15060 derived from oracle+cx_oracle://cgms12eu:***@eurdas.world
      smu_no=9050131, area=625000000, idstu=9000282 covering 50% of smu.
        Soil parameters {'SMLIM': 0.312, 'SMFCF': 0.312, 'SMW': 0.152, 'CRAIRC': 0.06,
                         'KSUB': 10.0, 'RDMSOL': 10.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.439}
      smu_no=9050131, area=625000000, idstu=9000283 covering 50% of smu.
        Soil parameters {'SMLIM': 0.28325, 'SMFCF': 0.28325, 'SMW': 0.12325, 'CRAIRC': 0.06,
                         'KSUB': 10.0, 'RDMSOL': 40.0, 'K0': 10.0, 'SOPE': 10.0, 'SM0': 0.42075}
    >>> for smu_no, area, idstu, percentage, soil_par in soildata:
    ...     print(smu_no, area, idstu, percentage)
    ...
    (9050131, 625000000, 9000282, 50)
    (9050131, 625000000, 9000283, 50)
    """

    def __init__(self, engine, idgrid, idcover=1000):
        # Initialise
        list.__init__(self)
        metadata = MetaData(engine)
        self.db_resource = str(engine)[7:-1]
        self.idgrid = int(idgrid)

        # Loop over all the mapping units
        SMUs = self._get_SMU_from_EMU(metadata, self.idgrid, idcover)
        for idsmu, area in SMUs:
            STUs = self._get_STU_from_SMU(metadata, idsmu)
            for idstu, percentage in STUs:
                soil_par = SoilDataProviderSingleLayer(engine, idstu)
                self.append((idsmu, area, idstu, percentage, soil_par))

    def _get_SMU_from_EMU(self, metadata, idgrid, idcover):
        """Retrieves the relevant SMU for given idgrid from table link_smu_grid_cover."""
        result = None
        t = Table("link_smu_grid_cover", metadata, autoload=True)
        sm = select([t.c.idsmu, t.c.area], and_(t.c.idgrid == self.idgrid, t.c.idcover == idcover))
        sc = sm.execute()
        result = sc.fetchall()
        sc.close()
        if result is None or len(result) == 0:
            msg = "No soil mapping units found for idgrid=%s" % idgrid
            raise exc.PCSEError(msg)
        return result

    def _get_STU_from_SMU(self, metadata, idsmu):
        """Retrieves the relevant soil typologic units for given idsmu from table link_smu_stu"""
        result = None
        t = Table("link_smu_stu", metadata, autoload=True)
        sc = select([t.c.idstu, t.c.percentage], t.c.idsmu == idsmu).execute()
        result = sc.fetchall()
        sc.close()
        if result is None:
            msg = "No soil typologic units found for idsmu=%i" % idsmu
            raise exc.PCSEError(msg)
        return result

    def __str__(self):
        result = "Soil data for grid_no=%i derived from %s\n" % (self.idgrid, self.db_resource)
        template = "  idsmu=%i, area=%.0f, idstu=%i covering %i%% of smu.\n    Soil parameters %s\n"
        for t in self:
            result += template % t
        return result


class CropDataProvider(dict):
    """Retrieves the crop parameters for the given idgrid, idcrop and year
    from the tables CROP_CALENDAR, CROP_PARAMETER_VALUE and VARIETY_PARAMETER_VALUE.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid: Integer grid ID, maps to the idgrid column in the table
    :param idcrop: Integer crop ID, maps to the idcrop column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """
    # Define single and tabular crop parameter values
    parameter_codes_single = wofost_parameters.WOFOST_parameter_codes_single
    parameter_codes_tabular = wofost_parameters.WOFOST_parameter_codes_tabular
    # Some parameters have to be converted from a single to a tabular form
    single2tabular = wofost_parameters.WOFOST_single2tabular
    # Default values for additional parameters not defined in CGMS
    parameters_additional = wofost_parameters.WOFOST_parameters_additional
    # Optional parameters, mainly dealing with vernalisation
    parameters_optional = wofost_parameters.WOFOST_optional_parameters

    def __init__(self, engine, idgrid, idcrop_parametrization):
        dict.__init__(self)
        self.idgrid = int(idgrid)
        self.idcrop_parametrization = int(idcrop_parametrization)
        self.crop_name = fetch_crop_name(engine, idcrop_parametrization)
        self.db_resource = str(engine)

        metadata = MetaData(engine)

        # Get crop variety from crop_spatializations
        t = Table('crop_spatializations', metadata, autoload=True)
        sm = select([t.c.idvariety],
                    and_(t.c.idgrid == self.idgrid,
                         t.c.idcrop_parametrization == self.idcrop_parametrization))
        sc = sm.execute()
        row = sc.fetchone()
        sc.close()
        if row is None:
            msg = ("No entry found in table crop_spatializations for idgrid=%s and "
                   "idcrop_parametrization=%s" % (self.idgrid, self.idcrop_parametrization))
            raise exc.PCSEError(msg)
        self.idvariety = int(row.idvariety)

        # get the parameters from the CGMS db
        self._fetch_crop_parameter_values(metadata, self.idcrop_parametrization)
        self._fetch_variety_parameter_values(metadata, self.idcrop_parametrization, self.idvariety)
        self.update(self.parameters_additional)

        # Finally add crop name
        self["CRPNAM"] = self.crop_name

    def _fetch_crop_parameter_values(self, metadata, idcrop_parametrization):
        """Derived the crop parameter values from the table crop_parametrization_parameter
        table for given idcrop and add directly to dict self[]..
        """
        t1 = Table("crop_parametrization_parameter", metadata, autoload=True)
        t2 = Table("global_crop_parameters", metadata, autoload=True)

        # Pull single value parameters from table crop_parametrization_parameter
        sm = select([t1.c.xvalue, t2.c.crop_parameter],
                    and_(t1.c.idcrop_parametrization == idcrop_parametrization,
                         t1.c.idcrop_parameter == t2.c.idcrop_parameter,
                         t2.c.multi == 'N', t2.c.idcategory == 1))
        sc = sm.execute()
        rows = sc.fetchall()
        sc.close()

        for row in rows:
            if row.crop_parameter not in self.single2tabular:
                self[row.crop_parameter] = float(row.xvalue)
            else:
                pvalue = float(row.xvalue)
                code, value = self._convert_single2tabular(row.crop_parameter, pvalue)
                self[code] = value

        # Check that we have had all the single and single2tabular parameters now
        for parameter_code in (self.parameter_codes_single + tuple(self.single2tabular.keys())):
            found = False
            if parameter_code not in self.single2tabular:
                if parameter_code in self:
                    found = True
            else:
                for key in self:
                    if key.startswith(parameter_code):
                        found = True
                        break
            if not found and parameter_code not in self.parameters_optional:
                msg = ("No parameter value found for idcrop_parametrization=%s, "
                       "parameter_code='%s'." % (self.idcrop_parametrization, parameter_code))
                raise exc.PCSEError(msg)

        # Pull tabular parameters from crop_parametrization_parameter
        for crop_parameter in self.parameter_codes_tabular:
            pattern = crop_parameter + r'%'
            sc = select([t1.c.xvalue, t1.c.yvalue, t2.c.crop_parameter],
                        and_(t1.c.idcrop_parametrization == idcrop_parametrization,
                             t1.c.idcrop_parameter == t2.c.idcrop_parameter,
                             t2.c.idcategory == 1, t2.c.multi == 'Y',
                             t2.c.crop_parameter.like(pattern)),
                        order_by=[t2.c.crop_parameter]).execute()
            rows = sc.fetchall()
            sc.close()
            if not rows and crop_parameter not in self.parameters_optional:
                msg = "No parameter value found for idcrop_parametrization=%s, crop_parameter='%s'."
                raise exc.PCSEError(msg % (self.idcrop_parametrization, crop_parameter))

            if len(rows) == 1:
                msg = ("Single parameter value found for idcrop_parametrization=%s, "
                       "crop_parameter='%s' while tabular parameter expected." %
                       (idcrop_parametrization, crop_parameter))
                raise exc.PCSEError(msg)
            values = []
            for row in rows:
                values.extend([float(row.xvalue), float(row.yvalue)])
            self[crop_parameter] = values

    def _fetch_variety_parameter_values(self, metadata, idcrop_parametrization, idvariety):
        """Derived the crop parameter values from the table crop_variety_parameters
        for given idvariety_on and add directly to dict self[].
        """
        t1 = Table("crop_variety_parameters", metadata, autoload=True)
        t2 = Table("global_crop_parameters", metadata, autoload=True)

        # Pull single value parameters from table crop_variety_parameter
        sm = select([t1.c.xvalue, t2.c.crop_parameter],
                    and_(t1.c.idvariety == idvariety,
                         t1.c.idcrop_parameter == t2.c.idcrop_parameter,
                         t2.c.multi == 'N', t2.c.idcategory == 1))
        sc = sm.execute()
        rows = sc.fetchall()
        sc.close()

        # Loop over the rows
        for row in rows:
            if row.crop_parameter in self.parameter_codes_single:
                self[row.crop_parameter] = float(row.xvalue)
            elif row.crop_parameter in self.single2tabular:
                pvalue = float(row.xvalue)
                code, value = self._convert_single2tabular(row.crop_parameter, pvalue)
                self[code] = value

        # Pull tabular parameters from crop_variety_parameter
        for crop_parameter in self.parameter_codes_tabular:
            pattern = crop_parameter + r'%'
            sc = select([t1.c.xvalue, t1.c.yvalue, t2.c.crop_parameter],
                        and_(t1.c.idvariety == idvariety,
                             t1.c.idcrop_parameter == t2.c.idcrop_parameter,
                             t2.c.multi == 'Y', t2.c.idcategory == 1,
                             t2.c.crop_parameter.like(pattern)),
                        order_by=[t2.c.crop_parameter]).execute()
            rows = sc.fetchall()
            sc.close()
            if not rows:
                continue
            if len(rows) == 1:
                msg = ("Single parameter value found for idcrop_parametrization=%s, "
                       "crop_parameter='%s' while tabular parameter expected."
                       % (idcrop_parametrization, crop_parameter))
                raise exc.PCSEError(msg)
            values = []
            for row in rows:
                values.extend([float(row.xvalue), float(row.yvalue)])
            self[crop_parameter] = values

    def _convert_single2tabular(self, crop_parameter, pvalue):
        """Converts the single parameter into a tabular parameter.
        """
        tabular_crop_parameter, template = self.single2tabular[crop_parameter]
        tabular_values = [pvalue if v is None else v for v in template]
        return tabular_crop_parameter, tabular_values

    def __str__(self):
        msg = ("Crop parameter values for idgrid=%s, idcrop_parametrization=%s (%s), "
               "idvariety=%s derived from %s\n" % (self.idgrid, self.idcrop_parametrization,
                                                      self.crop_name, self.idvariety,
                                                      self.db_resource))
        msg += str(self)
        return msg


class SiteDataProvider(dict):
    """Provides the site data from the tables INITIAL_SOIL_WATER and SITE.

    :param engine: SqlAlchemy engine object providing DB access
    :param idgrid:  Grid number (int)
    :param idcrop: Crop number (int)
    :param campaign_year: Campaign year (int)
    :param idstu: soil typologic unit number (int)

    Note that the parameter SSI, SSMAX, NOTIF, IFUNRN are not defined in the
    CGMS14 database and are set to zero.

    Moreover, the start date of the water balance is defined by the
    column POTENTIAL_WATER_STARTDATE. This value can be accessed as
    an attribute `start_date_waterbalance`.
    """
    _defaults = {"IFUNRN": 0,
                 "NOTINF": 0,
                 "SSMAX": 0.0,
                 "SSI": 0.0}

    def __init__(self, engine, idgrid, idcrop_parametrization, campaign_year, idstu):
        # Initialise
        dict.__init__(self)
        self.idgrid = int(idgrid)
        self.idcrop_parametrization = int(idcrop_parametrization)  # TODO: check!!
        self.campaign_year = int(campaign_year)
        self.idstu = int(idstu)

        # Start using the engine
        self.db_resource = str(engine)[7:-1]
        self.crop_name = fetch_crop_name(engine, self.idcrop_parametrization)
        metadata = MetaData(engine)

        # Now select from the appropriate table
        t = Table('soil_initial_water', metadata, autoload=True)
        sm = select([t],
                    and_(t.c.idgrid == self.idgrid, t.c.idcrop_parametrization == idcrop_parametrization,
                         t.c.year == self.campaign_year, t.c.idstu == self.idstu))
        sc = sm.execute()
        row = sc.fetchone()
        sc.close()
        if row is None:
            msg = ("Failed retrieving site data for grid_no=%s, idcrop_parametrization=%s, campaign_year=%s, "
                   "stu_no=%s" % (self.idgrid, self.idcrop_parametrization, self.campaign_year, self.idstu))
            raise exc.PCSEError(msg)

        # This is what we can learn from table soil_initial_water:
        self["WAV"] = float(row.rooting_depth_potential_water)

        # Start date water balance
        self.start_date_waterbalance = check_date(row.potential_water_startdate)

        # Set some missing parameters in CGMS14 to default values
        self.update(self._defaults)

    def __str__(self):
        result = ("Site parameter values for grid_no=%s, idcrop_parametrization=%s (%s), stu_no=%s, "
                  "campaign_year=%i derived from %s\n" % (self.idgrid, self.idcrop_parametrization,
                                                          self.crop_name, self.idstu, self.campaign_year,
                                                          self.db_resource))
        result += "    %s" % dict.__str__(self)
        return result
