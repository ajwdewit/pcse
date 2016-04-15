# -*- coding: utf-8 -*-
# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2016

"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers are compatible with a CGMS 14 database schema.
"""

import datetime

from sqlalchemy import MetaData, select, Table, and_
# from tabulate import tabulate
import numpy as np

from ...util import wind10to2, safe_float, check_date, reference_ET
from ... import exceptions as exc
from ...base_classes import WeatherDataContainer, WeatherDataProvider

def fetch_crop_name(engine, crop_no):
    """Retrieves the name of the crop from the CROP table for given crop_no.

    :param engine: SqlAlchemy engine object providing DB access
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """
    result = ""
    metadata = MetaData(engine)
    table_crop = Table("global_crops", metadata, autoload=True)
    sc = select([table_crop], table_crop.c.idcrop == crop_no).execute()
    row = sc.fetchone()
    sc.close()
    if row is None:
        msg = "Failed deriving crop name from view GLOBAL_CROPS for crop_no %s" % crop_no
        raise exc.PCSEError(msg)
    result = row.descriptor
    return result

class STU_Suitability(set):
    """Returns a set() of suitable STU's for given crop_no.

    :param engine: SqlAlchemy engine object providing DB access
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    """
    def __init__(self, engine, crop_no):
        self.crop_no = int(crop_no)

class WeatherObsGridDataProvider(WeatherDataProvider):
    """Retrieves meteodata from the WEATHER_OBS_GRID table in a CGMS14
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
    """
        # default values for the Angstrom parameters in the sunshine duration model
    angstA = 0.18
    angstB = 0.55
   
    def __init__(self, engine, grid_no, start_date=None, end_date=None, recalc_ET=False):        WeatherDataProvider.__init__(self)
        self.grid_no = grid_no
        self.recalc_ET = recalc_ET
        
    #---------------------------------------------------------------------------
    def _fetch_location_from_db(self, metadata):
        """Retrieves latitude, longitude, elevation from "grid" table and
        assigns them to self.latitude, self.longitude, self.elevation."""
       pass
   
    #---------------------------------------------------------------------------
    def _fetch_weather_from_db(self, metadata):
        """Retrieves the meteo data from table "grid_weather".
        """
        pass
    
    #---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, t):
        """Process record from grid_weather including unit conversion."""
        result = None
        
        return result
    
class TimerDataProvider(dict):
    """Class for providing timerdata from the CROP_CALENDAR table in a CGMS14 database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the GRID_NO column in the table
    :param crop_no: Integer crop ID, maps to the CROP_NO column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
        The campaign year usually refers to the year of the harvest. Thus for crops
        crossing calendar years, the start_date can be in the previous year.
    """
    def __init__(self, engine, grid_no, crop_no, campaign_year):   
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no) 
        self.campaign_year = int(campaign_year)
    
    def set_START_DATE(self, start_date):
        """Updates the value for START_DATE in TimerDataProvider
        """
        pass
    
    def __str__(self):
        result = ""
        
        return result

class SoilDataProviderSingleLayer(dict):
        """Class for providing soil data from the ROOING_DEPTH AND
    SOIL_PHYSICAL_GROUP tableS in a CGMS14 database. This
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
    
    def _get_from_STU(self, metadata, stu_no):
        """Retrieves the soil parameters for the given soil typologic unit
        (stu_no) from the tables SOIL_PHYSICAL_GROUP and ROOTING_DEPTH.
        """
        rd_class = 0
        soil_group_no = 0
        
        return (rd_class, soil_group_no)
    
    def _get_rooting_depth(self, metadata, rd_class):
        """Gets the rooting depth from the table ROOTING_DEPTH and
         stores into self[] directly under parameter name 'RDMSOL'.

        :param metadata: An SQLAlchemy Metadata object
        :param rd_class: The rooting depth class (integer)
        """
        pass
    
    def _get_soil_hydraulic_parameters(self, metadata, spg_no):
        """Retrieves the soil hydraulic parameters and stores into self[] directly.

        :param metadata: An SQLAlchemy Metadata object
        :param spg_no: the soil physical group number (integer)
        :return: None
        """
        pass
    
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
    
    def __init__(self, engine, grid_no):
        list.__init__(self)
        self.grid_no = int(grid_no)

    def _get_SMU_from_EMU(self, metadata, grid_no):
        """Retrieves the relevant SMU for given grid_no from
        table EMU.
        """
        result = None
        
        if result is None:
            msg = ("No soil mapping units (SMU) found for grid_no=%i in table EMU" % grid_no)
            raise exc.PCSEError(msg)
        return result
    
    def _get_STU_from_SMU(self, metadata, smu_no):
        """Retrieves the relevant STU for given SMU_NO from table
        SOIL_ASSOCIATION_COMPOSITION
        """
        result = None
        
        if result is None:
            msg = "No soil typologic units (STU) found for smu_no=%i" % smu_no
            raise exc.PCSEError(msg)
        return result
    
    def __str__(self):
        result = ""
        
        return result
    
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
        
    def _fetch_crop_parameter_values(self, metadata, crop_no):
        """Derived the crop parameter values from the CROP_PARAMETER_VALUE
        table for given crop_no and add directly to dict self[]..
        """
        pass
    
    def _fetch_variety_parameter_values(self, metadata, crop_no, variety_no):
        """Derived the crop parameter values from the VARIETY_PARAMETER_VALUE
        table for given crop_no & variety_on and add directly to dict self[].
        """
        pass
    
    def _convert_single2tabular(self, parameter_code, pvalue):
        """Converts the single parameter into a tabular parameter.
        """
        tabular_parameter_code = ""
        tabular_values = []
        
        return tabular_parameter_code, tabular_values
    
    def __str__(self):
        result = ""
        
        return result
	
class SiteDataProvider(dict):
    """Provides the site data from the tables INITIAL_SOIL_WATER and SITE.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no:  Grid number (int)
    :param crop_no: Crop number (int)
    :param campaign_year: Campaign year (int)
    :param stu_no: soil typologic unit number (int)

    Note that the parameter SSI (Initial surface storage) is
    set to zero

    """
    
    def __init__(self, engine, grid_no, crop_no, campaign_year, stu_no):
        dict.__init__(self)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.stu_no = int(stu_no)

    
    def __str__(self):
        result = ""
        
        return result