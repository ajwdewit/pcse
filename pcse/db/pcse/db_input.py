# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Routines for retrieving data from the PCSE demo database.
 
Implements the following functions:
    - fetch_cropdata()
    - fetch_sitedata()
    - fetch_soildata()
    - fetch_timerdata()


Implements the following classes:
    - MeteoFetcher

Implements the following exceptions:
    - MeteodataError
    - CropdataError
    - SitedataError
    - SoildataError
    - TimerdataError
"""
import sys, os
import datetime as dt
import logging
import copy
import array
from math import log10

from sqlalchemy import create_engine, MetaData, select, Table, and_, join
import yaml

from ...util import wind10to2, Afgen, check_date
from ... import exceptions as exc
from ...base_classes import WeatherDataContainer, WeatherDataProvider


#-------------------------------------------------------------------------------
class MeteodataError(exc.PCSEError):
    """Exception for EnsembleMeteoCl."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
         return repr(self.value)

#----------------------------------------------------------------------------
class CropdataError(exc.PCSEError):
    def __init__(self, value):
        self.value = value
    def __str__(self):
         return repr(self.value)

#-------------------------------------------------------------------------------
class SitedataError(exc.PCSEError):
    def __init__(self, value):
        self.value = value
    def __str__(self):
         return repr(self.value)

#-------------------------------------------------------------------------------
class SoildataError(exc.PCSEError):
    def __init__(self, grid, msg):
        value = "Failed to select soil type for grid %s: %s" % (grid, msg)
        self.value = value
    def __str__(self):
         return repr(self.value)

#-------------------------------------------------------------------------------
def fetch_cropdata(metadata, grid, year, crop):
    """Retrieve crop parameter values for given grid, year, crop from DB.
    
    Parameter values are pulled from tables 'crop_parameter_value'
    and 'variety_parameter_value'. Metadata is an SQLAlchemy metadata object.
    
    Returns a dictionary with WOFOST crop parameter name/value pairs.
    
    Note that the parameter names are defined in the function itself in order
    to distinguish scalar and table parameters. These definitions will need to
    be extended when additional parameters are to be retrieved.
    """
    
    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger('PCSE.db_util')

    # Create initial dictionary 
    cropdata={}    

    # Get crop name from crop table;
    table_cc = Table('crop', metadata, autoload=True)
    s = select([table_cc], table_cc.c.crop_no==crop)
    row = s.execute().fetchone()
    cropdata["CRPNAM"] = row.crop_name

    # Get crop variety from crop_calendar;
    table_cc = Table('crop_calendar', metadata, autoload=True)
    s = select([table_cc],
               and_(table_cc.c.grid_no==grid,
                    table_cc.c.crop_no==crop,
                    table_cc.c.year==year))
    r = s.execute()
    rows = r.fetchall()
    r.close()
    if (rows is not None) and (len(rows) == 1):
        variety = rows[0].variety_no
    else:
        logger.error(("None or multiple crop definitions found for grid: %7i and crop " + \
                            "no: %3i and year: %5i") % (grid, crop, year))
        raise CropdataError(("No crop definition found for grid: %7i and crop " + \
                            "no: %3i and year: %5i") % (grid, crop, year))
    
    # Define crop parameter values
    parameter_codes_sngl = ("CFET", "CVL", "CVO", "CVR", "CVS", "DEPNR", "DLC", 
                            "DLO", "DVSEND", "EFF", "IAIRDU", "IDSL", "KDIF", 
                            "LAIEM", "PERDL", "Q10", "RDI", "RDMCR", "RGRLAI", 
                            "RML", "RMO", "RMR", "RMS", "RRI", "SPA", "SPAN", "SSA", 
                            "TBASE", "TBASEM", "TDWI", "TEFFMX", "TSUM1", "TSUM2", 
                            "TSUMEM", "IOX")
    parameter_codes_mltp = ("AMAXTB", "DTSMTB", "FLTB", "FOTB", "FRTB", "FSTB", 
                            "RDRRTB", "RDRSTB", "RFSETB", "SLATB", "TMNFTB", 
                            "TMPFTB")
    
    # Pull single value parameters from CROP_PARAMETER_VALUE first
    table_crop_pv = Table('crop_parameter_value', metadata, autoload=True)
    for paramcode in parameter_codes_sngl:
        s = select([table_crop_pv], and_(table_crop_pv.c.crop_no==crop,
                                         table_crop_pv.c.parameter_code==paramcode))
        r = s.execute()
        rows = r.fetchall()
        if (len(rows) == 0):
            logger.error(("No crop parameter value found for:" + \
                          "crop %s, parameter: %s") % (crop, paramcode))
            raise CropdataError(("No crop parameter value found for:" + \
                                "crop %s, parameter: %s") % (crop, paramcode))
        elif (len(rows) == 1):
            cropdata[paramcode] = float(rows[0].parameter_xvalue)
        else:
            logger.error(("Multiple crop parameter values found for:" + \
                          "crop %s, parameter: %s") % (crop, paramcode))
            raise CropdataError(("Multiple crop parameter values found for:" + \
                                "crop %s, parameter: %s") % (crop, paramcode))
    logger.debug("Succesfully retrieved single value parameters "+\
                 "from CROP_PARAMETER_VALUE TABLE")
    
    # Pull array parameter values from CROP_PARAMETER_VALUE
    # note the change in the mask value and the use of "LIKE" in the SQL query
    for paramcode in parameter_codes_mltp:
        pattern = paramcode + r'%'
        s = select([table_crop_pv],
                   and_(table_crop_pv.c.crop_no == crop,
                        table_crop_pv.c.parameter_code.like(pattern)),
                   order_by=[table_crop_pv.c.parameter_code])
        r = s.execute()
        rows = r.fetchall()
        c = len(rows)
        if (c == 0):
            logger.error(("No crop parameter value found for:" + \
                          "crop %s, parameter: %s") % (crop, paramcode))
            raise CropdataError(("No crop parameter value found for:" + \
                                "crop %s, parameter: %s") % (crop, paramcode))
        elif (c == 1):
            logger.error(("Single crop parameter value found for:" + \
                          "crop %s, parameter: %s") % (crop, paramcode))
            raise CropdataError(("Single crop parameter value found for:" + \
                                "crop %s, parameter: %s") % (crop, paramcode))
        else:
            value = array.array('d', [0.]*(c*2))
            for i in range(0,c):
                value[i*2] = float(rows[i].parameter_xvalue)
                value[(i*2)+1]= float(rows[i].parameter_yvalue)
            cropdata[paramcode] = value
    logger.debug("Succesfully retrieved array parameters "+\
                 "from CROP_PARAMETER_VALUE TABLE")

    # Pull same parameter values from VARIETY_PARAMETER_VALUES
    # if they are defined for that variety.
    # Pull single value parameters first
    table_var_pv = Table('variety_parameter_value', metadata, autoload=True)
    for paramcode in parameter_codes_sngl:
        s = select([table_var_pv],
                   and_(table_var_pv.c.variety_no==variety,
                        table_var_pv.c.crop_no==crop,
                        table_var_pv.c.parameter_code==paramcode))
        r = s.execute()
        rows = r.fetchall()
        c = len(rows)
        if (c == 0):
            pass
        elif (c == 1):
            cropdata[paramcode] = float(rows[0].parameter_xvalue)
        else:
            errstr = "Multiple values found for: crop: %s, variety: %s, " + \
                     "parameter: %s which is supposed to be a single value."
            logger.error(errstr % (crop, paramcode))
            raise CropdataError(errstr % (crop, paramcode))
    logger.debug("Succesfully retrieved single value parameters "+\
                 "from VARIETY_PARAMETER_VALUE TABLE")
            
    # pull array value parameters - note the change in the mask value and
    # the use of "LIKE" in the SQL query
    for paramcode in parameter_codes_mltp:
        pattern = paramcode + r'%'
        s = select([table_var_pv],
                   and_(table_var_pv.c.crop_no==crop,
                        table_var_pv.c.parameter_code.like(pattern),
                        table_var_pv.c.variety_no==variety),
                   order_by=[table_var_pv.c.parameter_code])
        r = s.execute()
        rows = r.fetchall();
        c = len(rows)
        if (c == 0):
            pass
        elif (c == 1):
            errstr = "Single value found for: crop: %s, variety: %s, " + \
                     "parameter: %s which is supposed to be a single value."
            logger.error(errstr % (crop, paramcode))
            raise CropdataError(errstr % (crop, paramcode))
        else:
            value = array.array('d', [0.]*(c*2))
            for i in range(0,c):
                value[i*2] = float(rows[i].parameter_xvalue)
                value[(i*2)+1]= float(rows[i].parameter_yvalue)
            cropdata[paramcode] = value
    logger.debug("Succesfully retrieved array parameters "+\
                 "from VARIETY_PARAMETER_VALUE TABLE")

    # Make some specific changes for FORTRAN wofost with regard to variables
    # SSA, KDIF and EFF. This is needed because the FORTRAN code expects a
    # parameter array, while these parameters have been defined in CGMS as
    # single values. DVSI does not exist in CGMS, therefore set DVSI to zero.
    
    # SSA convert to SSATB:
    SSA = cropdata["SSA"]
    SSATB = array.array('d',[0.]*4)
    SSATB[1] = SSA
    SSATB[2] = 2.0
    SSATB[3] = SSA
    cropdata.update({"SSATB":SSATB})
    # KDIF convert to KDIFTB:
    KDIF = cropdata["KDIF"]
    KDIFTB = array.array('d',[0.]*4)
    KDIFTB[1] = KDIF
    KDIFTB[2] = 2.0
    KDIFTB[3] = KDIF
    cropdata.update({"KDIFTB":KDIFTB})
    # EFF convert to EFFTB
    EFF = cropdata["EFF"]
    EFFTB = array.array('d',[0.]*4)
    EFFTB[1] = EFF
    EFFTB[2] = 40.0
    EFFTB[3] = EFF
    cropdata.update({"EFFTB":EFFTB})
    # DVSI set to 0
    cropdata.update({"DVSI":0})
    
    logger.info("Succesfully retrieved crop parameter values from database")
    return cropdata

#-------------------------------------------------------------------------------
def fetch_soilparams(metadata, grid, soilgroup):
    """Retrieve soil parameters for given soilgroup from table SOIL_PHYSICAL_GROUP.
    
    Returns a dictionary with WOFOST soil parameter name/value pairs.
    
    Note that the grid parameter is only passed because it is needed when
    raising a SoildataError
    """
    
    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger('PCSE41.db_input') #Logger()

    # Define soil physical variable parameter codes
    # defined as (code_parname, db_parname)
    soil_parameters = [("CRAIRC", "CRITICAL_AIR_CONTENT"),
                       ("K0", "HYDR_CONDUCT_SATUR"), 
                       ("SOPE", "MAX_PERCOL_ROOT_ZONE"),
                       ("KSUB", "MAX_PERCOL_SUBSOIL"),
                       ("SMFCF", "SOIL_MOISTURE_CONTENT_FC"),
                       ("SM0", "SOIL_MOISTURE_CONTENT_SAT"),
                       ("SMW", "SOIL_MOISTURE_CONTENT_WP")]
    # Table soil parameters can be mapped directly to parameter names
    # in the WOFOST code.
    soil_parameters_mltp = ["CONTAB", "SMTAB"]

    # Select soil properties from table.
    soilparams = {}
    table_soil_pg = Table('soil_physical_group',metadata, autoload=True)
    for (code_parname, db_parname) in soil_parameters:
        r = select([table_soil_pg], 
                   and_(table_soil_pg.c.soil_group_no==soilgroup,
                        table_soil_pg.c.parameter_code==db_parname)).execute()
        rows = r.fetchall()
        if len(rows) == 0:
            msg = "missing parameter: %s" % db_parname
            raise SoildataError(grid, msg)
        soilparams[code_parname] = float(rows[0].parameter_xvalue)
            
        # as params or calculate???
        #     SMFCF(il)  = AFGEN(Soil(is)%SMfromPF, Soil(is)%ilSM, PFFC)       
        #     SMW(il)    = AFGEN(Soil(is)%SMfromPF, Soil(is)%ilSM, PFWP)
        #     SM0(il)    = Soil(is)%SMsat 
        #     CRAIRC(il) = Soil(is)%CritAir

    # Select soil physical array properties with like (table parameters)
    for soil_par in soil_parameters_mltp:
        pattern = soil_par + r'%'
        r = select([table_soil_pg], 
                   and_(table_soil_pg.c.soil_group_no==soilgroup,
                        table_soil_pg.c.parameter_code.like(pattern)),
                   order_by=[table_soil_pg.c.parameter_code]).execute()
        rows = r.fetchall()
        nrows = len(rows)
        if (nrows < 3):
            msg = ("Contains less than 3 (x,y) pairs for:" + 
                   "soilgroup %s, parameter: %s" % (soilgroup, soil_par))
            raise SoildataError(grid, msg)

        value =  array.array('d', [0.]*(nrows*2))
        for i, row in enumerate(rows):
            value[i*2] = float(row.parameter_xvalue)
            value[(i*2)+1]= float(row.parameter_yvalue)
        soilparams[soil_par] = Afgen(value)     # @ToDo: bij definitie AfgenTrait() zie assimilation.py

        if (soil_par == "SMTAB"):
            smtab = value
            ilsm = nrows
            # saturated water content; porosity
            # soilparams["SM0"] = soilparams["SMTAB"](-1.0);
            # SMsat, soil porosity, SM at Pf = -1.0 (0.0 ?)
            # formation of PFTAB (PFfromSM), the inverse of SMTAB (SMfromPF)
            value =  array.array('d', [0.]*(nrows*2))
            for i, row in enumerate(reversed(rows)):
                value[i*2]    = float(row.parameter_yvalue)
                value[(i*2)+1]= float(row.parameter_xvalue)
            soilparams["PFTAB"] = Afgen(value)  # @ToDo: bij definitie AfgenTrait() zie assimilation.py

            # Gauss points and weights
            NGAU = 5
            PGAU = array([0.0469100770, 0.2307653449, 0.5, 0.7692346551, \
                             0.9530899230])         # sum(PGAU) = 2.5?
            WGAU = array([0.1184634425, 0.2393143352, 0.2844444444, 0.2393143352, \
                             0.1184634425])         # sum(WGAU) = 0.9999999998 ~ 1
            ELOG10 = 2.302585092994
            
            # Calculation of Matric Flux Potential (MFP) for the soil moisture
            # content values SM given in interpolation function SMTAB. The MFP is
            # defined as the integral of the conductivity K(teta) from -Infinity to
            # a certain teta (teta is the soil moisture content SM).
            # The MPF is calculated as an integral of (K(teta) * 10^(pF) * eLog10)
            # over the pF range considered.
            value = array.array('d', [0.]*(ilsm*2))
            
            # set zero potential at highest PF
            value[ilsm-2] = smtab[ilsm-2]
            value[ilsm-1] = 0.0
            
            for i in range(ilsm-3, -1, -2):
                # copy PF-value as x-value
                value[i-1] = smtab[i-1]
                #get integral over PF range
                ADD = 0.0
                DeltaPF = smtab[i+1] - smtab[i-1]
                for j in range(0, NGAU):
                    PFg  = smtab[i-1] + PGAU[j] * DeltaPF
                    CON  = pow(10, soilparams["CONTAB"](PFg))
                    ADD += CON * 10**PFg * ELOG10 * WGAU[j]
                value[i]= value[i+2] + ADD * DeltaPF
            soilparams["MFPTAB"] = Afgen(value)  # @ToDo: bij definitie AfgenTrait() zie assimilation.py
            
            EquilTableLEN = 30
            value         = array.array('d', [0.]*EquilTableLEN) # WaterFromHeight
            inversevalue  = array.array('d', [0.]*EquilTableLEN) # HeightFromAir
            # Calculate soil air volume above watertable at equilibrium:
            # (this property is needed only for the subsoil type for groundwater run,
            # but is calculated for all types and situations)
            # ---
            # Table WaterFromHeight is calculated, containing the cumulative amount
            # of water as a function of height above groundwater under equilibrium
            # conditions. The table is calculated for the following values of height
            # (= matric head) : 0,2,4,8,16,32,64,128,256,.....16384.
            # method is 5 point gauss-lagrange integration of SM on each interval.
            # ---
            # Table HeightFromAir is found by inverting WaterFromHeight (and
            # subtracting water from pore space)
            MH0 = 0.0
            MH1 = 2.0
            value[0] = 0.0
            value[1] = 0.0
            inversevalue[0] = 0.0
            inversevalue[1] = 0.0
            
            for i in range (2, EquilTableLEN, 2):
                value[i]   = MH1
                value[i+1] = value[i-1]
                for i3 in range(0, NGAU):
                    PFg = log10(MH0 + PGAU[i3] * (MH1-MH0))
                    value[i+1] += WGAU[i3] * (MH1-MH0) * soilparams["SMTAB"](PFg)
                
                # and the inverse
                inversevalue[i]   = value[i]*soilparams["SM0"] - value[i+1]
                inversevalue[i+1] = value[i]
                
                # new interval
                MH0 = MH1
                MH1 = 2*MH1
            soilparams["WaterFromHeight"] = Afgen(value)        # @ToDo: bij definitie AfgenTrait() zie assimilation.py
            soilparams["HeightFromAir"] = Afgen(inversevalue)   # @ToDo: bij definitie AfgenTrait() zie assimilation.py

    logger.info("Succesfully retrieved soil parameter values from database")
    return soilparams

def fetch_soildata_layered(metadata, grid):
    """Retrieve soil parameters for a layered soil for given grid from DB.
    
    Retrieves the soil type from the table SOIL_TYPE and soil
    parameters from tables SOIL_LAYERS and SOIL_PHYSICAL_GROUP.
    
    Returns a dictionary with WOFOST soil parameter name/value pairs.
    """
    
    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger('PCSE.db_util')
    
    soildata = dict()
    table_soiltype = Table('soil_type', metadata, autoload=True)
    r = select([table_soiltype.c.grid_no,
                table_soiltype.c.soil_type_no],
                table_soiltype.c.grid_no==grid).execute()
    rows = r.fetchall()
    if len(rows) != 1:
        msg = "None or multiple records found in table SOIL_TYPE"
        raise SoildataError(grid, msg)
    soil_type_no = (rows[0]).soil_type_no

    # Select soil layer(s) from the table SOIL_LAYERS
    table_soillayer = Table('soil_layers', metadata, autoload=True)
    r = select([table_soillayer.c.thickness,
                table_soillayer.c.soil_group_no],
                table_soillayer.c.soil_type_no==soil_type_no,
                order_by=[table_soillayer.c.layer_no]).execute()
    rows = r.fetchall()
    if len(rows) == 0:
        msg = "No soil layers found for soil_type_no %i." % soil_type_no
        raise SoildataError(grid, msg)
    soildata["NSL"] = len(rows)

    # layered soil, get the layers
    soil_layers = []
    lbsl = 0.0
    for (thickness, soilgroup) in rows:
        lbsl += float(thickness)
        soil_params = fetch_soilparams(metadata, grid, soilgroup)
        layer_info = {"TSL":float(thickness),
                      "SOIL_GROUP_NO":soilgroup,
                      "LBSL":lbsl,
                      "SOILTYPE":soil_params}
        soil_layers.append(layer_info)
        
    soildata["SOIL_LAYERS"]= soil_layers
    
    #print "Succesfully retrieved soil layer parameter values for layered soil from database"
    logger.info("Succesfully retrieved soil layer parameter values for layered soil from database")

    return soildata

#-------------------------------------------------------------------------------
def fetch_soildata(metadata, grid):
    """Retrieve soil parameters for given grid from DB for a 1-layer soil.
    
    Retrieves soil_type_no from the table SOIL_TYPE and associated soil layers
    and soil physical data from tables SOIL_LAYERS and SOIL_PHYSICAL_GROUP. 
    
    Returns a dictionary with WOFOST soil parameter name/value pairs.
    """
    
    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger('PCSE.db_util')
    
    soildata = {}
    # Select soil from the table SOIL_TYPE
    table_soiltype = Table('soil_type', metadata, autoload=True)
    r = select([table_soiltype.c.grid_no,
                table_soiltype.c.soil_type_no],
                table_soiltype.c.grid_no==grid).execute()
    row = r.fetchone()
    r.close()
    if row is None:
        raise SoildataError(grid, "No record found!")
    soil_type_no = row.soil_type_no
    
    # Derive layers for given soil_type_no. This should return only one
    # layer, otherwise raise an error.
    table_soil_layers = Table('soil_layers',metadata, autoload=True)
    cursor = select([table_soil_layers.c.thickness,
                     table_soil_layers.c.soil_group_no],
                     table_soil_layers.c.soil_type_no==soil_type_no,
                     order_by=[table_soil_layers.c.layer_no]).execute()
    rows = cursor.fetchall()
    cursor.close()
    if len(rows) == 0:
        msg = "No record found."
        raise SoildataError(grid, msg)
    elif len(rows) > 1:
        msg = ("Number of soil layers > 1. Not possible for unlayered " +
               "waterbalance module. Use 'fetch_soiltype_multilayer'") 
        raise SoildataError(grid, msg)
    else:
        firstrow = rows[0]
    soildata["RDMSOL"] = float(firstrow[0])
    soil_group_no = firstrow[1]
    
    # Retrieve soil physical properties for given layer for given soil
    # parameter codes: (wofost_parname, database_name)
    soil_parameters = [("CRAIRC", "CRITICAL_AIR_CONTENT"),
                       ("K0", "HYDR_CONDUCT_SATUR"),
                       ("SOPE", "MAX_PERCOL_ROOT_ZONE"),
                       ("KSUB", "MAX_PERCOL_SUBSOIL"),
                       ("SMFCF", "SOIL_MOISTURE_CONTENT_FC"),
                       ("SM0", "SOIL_MOISTURE_CONTENT_SAT"),
                       ("SMW", "SOIL_MOISTURE_CONTENT_WP")]
    table_soil_pg = Table('soil_physical_group',metadata, autoload=True)
    for (wofost_soil_par, db_soil_par) in soil_parameters:
        r = select([table_soil_pg], 
                   and_(table_soil_pg.c.soil_group_no==soil_group_no,
                        table_soil_pg.c.parameter_code==db_soil_par)).execute()
        row = r.fetchone()
        if row is None:
            msg = "Parameter %s not found" % db_soil_par
            raise SoildataError(grid, msg)
        soildata[wofost_soil_par] = float(row.parameter_xvalue)

    logger.info("Succesfully retrieved soil parameter values from database")
    return soildata


class AgroManagementDataProvider(list):
    """Class for providing agromanagement data from the CROP_CALENDAR table in a PCSE database.

    :param engine: SqlAlchemy engine object providing DB access
    :param grid_no: Integer grid ID, maps to the grid_no column in the table
    :param crop_no: Integer id of crop, maps to the crop_no column in the table
    :param campaign_year: Integer campaign year, maps to the YEAR column in the table.
           The campaign year refers to the year of the crop start. Thus for crops
           crossing calendar years, the start_date can be in the previous year as the
           harvest.
    
    Note that this AgroManagementDataProvider is only used for the internal PCSE database
    and not to be used for CGMS databases.
    """
    agro_management_template = """
          - {campaign_start_date}:
                CropCalendar:
                    crop_id: '{crop_name}'
                    crop_start_date: {crop_start_date}
                    crop_start_type: {crop_start_type}
                    crop_end_date: {crop_end_date}
                    crop_end_type: {crop_end_type}
                    max_duration: {duration}
                TimedEvents: null
                StateEvents: null
        """

    def __init__(self, engine, grid_no, crop_no, campaign_year):
        list.__init__(self)
        self.grid_no = int(grid_no)
        self.crop_no = int(crop_no)
        self.campaign_year = int(campaign_year)
        self.crop_name = "not defined"

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

        # Determine the start date.
        self.crop_start_date = check_date(row.crop_start_date)
        self.campaign_start_date = self.crop_start_date

        # Determine the start date/type. Only sowing|emergence is accepted by PCSE/WOFOST
        self.crop_start_type = str(row.crop_start_type).strip()
        if self.crop_start_type not in ["sowing","emergence"]:
            msg = "Unrecognized crop start type: %s" % self.crop_start_type
            raise exc.PCSEError(msg)

        # Determine maximum duration of the crop
        self.max_duration = int(row.max_duration)

        # Determine crop end date/type and the end of the campaign
        self.crop_end_type = str(row.crop_end_type).strip().lower()
        if self.crop_end_type not in ["harvest", "earliest", "maturity"]:
            msg = ("Unrecognized option for END_TYPE in table "
                   "CROP_CALENDAR: %s" % row.end_type)
            raise exc.PCSEError(msg)

        if self.crop_end_type == "maturity":
            self.crop_end_date = "null"
        else:
            self.crop_end_date = check_date(row.crop_end_date)

        input = self._build_yaml_agromanagement()
        self._parse_yaml(input)

    def _build_yaml_agromanagement(self):
        """Builds the YAML agromanagent string"""

        input = self.agro_management_template.format(campaign_start_date=self.campaign_start_date,
                                                     crop_name=self.crop_no,
                                                     crop_start_date=self.crop_start_date,
                                                     crop_start_type=self.crop_start_type,
                                                     crop_end_date=self.crop_end_date,
                                                     crop_end_type=self.crop_end_type,
                                                     duration=self.max_duration
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


#-------------------------------------------------------------------------------
def fetch_sitedata(metadata, grid, year):
    """Retrieve site data from DB for given grid, year.
    
    Pulls sitedata from the PCSE database 'SITE' table,
    
    Returns a dictionary with site parameter name/value pairs.
    """

    # Define a logger for the PCSE db_util routines
    logger = logging.getLogger('PCSE.db_util')

    try:
        #Get all settings from table 'SITE'
        table_site = Table('site', metadata, autoload=True)
        r = select([table_site],
                   and_(table_site.c.grid_no==grid,
                        table_site.c.year==year)
                   ).execute()
        row = r.fetchone()
        r.close()
        if row is not None:
            sitedata = {}
            sitedata['IFUNRN'] = float(row.ifunrn)
            sitedata['SSMAX']  = float(row.max_surface_storage)
            sitedata['NOTINF'] = float(row.not_infiltrating_fraction)
            sitedata['SSI']    = float(row.initial_surface_storage)
            sitedata['WAV']    = float(row.inital_water_availability)
            sitedata['SMLIM']  = float(row.smlim)
        else:
            raise RuntimeError("No rows found")
    except Exception as e:
        errstr = "Failed to get site data for year %i and grid %i: " + str(e)
        logger.error(errstr % (year, grid))
        raise SitedataError(errstr % (year, grid))

    logger.info("Successfully retrieved site variables from database")
    return sitedata


#----------------------------------------------------------------------------
class EnsembleGridWeatherDataProvider(WeatherDataProvider):
    """Retrieves ensemble meteodata from database.
    
    :param metadata: SqlAlchemy metadata object providing DB access
    :param grid_no:  CGMS Grid ID
    :param startdate: Retrieve meteo data starting with startdate
        (datetime.date object)
    :param enddate: Retrieve meteo data up to and including enddate
        (datetime.date object)

    This class joins the tables 'ensemble_grid_weather' and 'grid_weather' to
    retrieve a complete set of weather variables for each ensemble member.
    
    Currently this WeatherDataProvider is configured to read ensemble
    rainfall data from the table 'ensemble_grid_weather' and read the other
    deterministic elements from the table 'grid_weather'. If other ensemble
    meteorologic elements are needed then the function
    `_fetch_ensemble_weather_from_db` needs to be adapted.
    
    Note that all meteodata is first retrieved from the DB and stored
    internally. Therefore, no DB connections are stored within the class
    instance. This makes that class instances can be pickled.    
    """
    supports_ensembles = True
    
    def __init__(self, metadata, grid_no, startdate, enddate):

        WeatherDataProvider.__init__(self)

        self.grid_no = grid_no
        self.startdate = startdate
        self.enddate = enddate
        self.timeinterval = (enddate - startdate).days + 1
        
        # Get location info and return a template WeatherDataContainer
        wdc = self._fetch_location_from_db(metadata)

        # Retrieve ensemble weather data
        self._fetch_ensemble_weather_from_db(metadata, wdc)
            
    #---------------------------------------------------------------------------
    def _fetch_location_from_db(self, metadata):
        """Retrieves latitude, longitude, elevation from 'grid' table and
        returns a template WeatherDataContainer with these values.
        """

        # Pull Latitude value for grid nr from database

        try:
            table_grid = Table('grid', metadata, autoload=True)
            r = select([table_grid.c.latitude, table_grid.c.longitude,
                        table_grid.c.altitude],
                       table_grid.c.grid_no==self.grid_no).execute()
            row = r.fetchone()
            r.close()
            if row is None:
                raise Exception
        except Exception as exc:
            msg = "Failed deriving location info for grid %s" % self.grid_no
            raise MeteodataError(msg)

        wdc = WeatherDataContainer(LAT=row.latitude, LON=row.longitude,
                                   ELEV=row.altitude)

        msg = "Succesfully retrieved location information from 'grid' table "+\
              "for grid %s"
        self.logger.info(msg % self.grid_no)
        
        return wdc
    
    #---------------------------------------------------------------------------
    def _fetch_ensemble_weather_from_db(self, metadata, wdc):
        """Retrieves the ensemble meteo data by combining tables
        'ensemble_grid_weather' and 'grid_weather'.
        """
        
        # Load table definitions
        gw = Table('grid_weather', metadata, autoload=True)
        egw = Table('ensemble_grid_weather', metadata, autoload=True)
        
        # Define which columns to take from which tables
        columns = [egw.c.grid_no, egw.c.day, egw.c.member_id,
                   egw.c.rainfall, gw.c.maximum_temperature,
                   gw.c.minimum_temperature, gw.c.calculated_radiation,
                   gw.c.vapour_pressure, gw.c.windspeed, gw.c.et0,
                   gw.c.es0, gw.c.e0]
        whereclause = and_(egw.c.grid_no==gw.c.grid_no,
                           egw.c.day==gw.c.day,
                           egw.c.grid_no==self.grid_no,
                           egw.c.day>=self.startdate,
                           egw.c.day<=self.enddate)

        # construct query and fetch all rows
        sjoined = select(columns, whereclause)
        rows = sjoined.execute().fetchall()
        
        # Process rows into WeatherDataContainers
        distinct_member_ids = []
        for row in rows:
            keydate = self.check_keydate(row.day)
            member_id = row.member_id
            twdc = self._make_WeatherDataContainer(row, wdc)
            self._store_WeatherDataContainer(twdc, keydate, member_id)
            
            if member_id not in distinct_member_ids:
                distinct_member_ids.append(member_id)
            
        # Check if number of rows matches the number of distinct members times
        # the length of the time series
        nmember_ids = len(distinct_member_ids)
        expected_rows = self.timeinterval * nmember_ids
        nrows = len(rows)
        if nrows < expected_rows:
            msg =  "Only %i records selected while %i expected "+\
                       "for grid %i (%i days X %i members)"
            self.logger.warn(msg % (nrows, expected_rows, self.grid_no,
                                    self.timeinterval, nmember_ids))
        else:
            msg = "Retrieved ensemble weather data for grid %s between "+\
                  "%s and %s for %i members"
            self.logger.info(msg % (self.grid_no, self.startdate,
                                    self.enddate, nmember_ids))

    #---------------------------------------------------------------------------
    def _make_WeatherDataContainer(self, row, wdc):
        """Process record from grid_weather including unit conversion."""

        twdc = copy.deepcopy(wdc)
        twdc.DAY = row.day
        twdc.add_variable("TMAX", float(row.maximum_temperature),"Celsius")
        twdc.add_variable("TMIN", float(row.minimum_temperature),"Celsius")
        twdc.add_variable("VAP",  float(row.vapour_pressure),"hPa")
        twdc.add_variable("WIND", wind10to2(float(row.windspeed)),"m/sec")
        twdc.add_variable("RAIN", float(row.rainfall)/10.,"cm/day")
        twdc.add_variable("E0",  float(row.e0)/10.,"cm/day")
        twdc.add_variable("ES0", float(row.es0)/10.,"cm/day")
        twdc.add_variable("ET0", float(row.et0)/10.,"cm/day")
        twdc.add_variable("IRRAD", float(row.calculated_radiation)*1000.,"J/m2/day")
        
        return twdc

