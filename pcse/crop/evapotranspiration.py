# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
from math import exp

import array
import numpy as np

from ..traitlets import Float, Int, Instance, Bool, Instance
from ..decorators import prepare_rates, prepare_states
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
                         SimulationObject
from ..util import limit, merge_dict, AfgenTrait


def SWEAF(ET0, DEPNR):
    """Calculates the Soil Water Easily Available Fraction (SWEAF).

    :param ET0: The evapotranpiration from a reference crop.
    :param DEPNR: The crop dependency number.
    
    The fraction of easily available soil water between field capacity and
    wilting point is a function of the potential evapotranspiration rate
    (for a closed canopy) in cm/day, ET0, and the crop group number, DEPNR
    (from 1 (=drought-sensitive) to 5 (=drought-resistent)). The function
    SWEAF describes this relationship given in tabular form by Doorenbos &
    Kassam (1979) and by Van Keulen & Wolf (1986; p.108, table 20)
    http://edepot.wur.nl/168025.
    """
    A = 0.76
    B = 1.5
    # curve for CGNR 5, and other curves at fixed distance below it
    sweaf = 1./(A+B*ET0) - (5.-DEPNR)*0.10

    # Correction for lower curves (CGNR less than 3)
    if (DEPNR < 3.):
        sweaf += (ET0-0.6)/(DEPNR*(DEPNR+3.))

    return limit(0.10, 0.95, sweaf)


class EvapotranspirationWrapper(SimulationObject):
    """This selects the appropriate evapotranspiration module
    depending on the use of a waterbalance with a layered or a non-layered
    soil and whether the impact of CO2 should be taken into account
    """
    etmodule = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: parameter values
        """

        if "soil_profile" in parvalues:
            self.etmodule = EvapotranspirationCO2Layered(day, kiosk, parvalues)
        elif "CO2" in parvalues and "CO2TRATB" in parvalues:
            self.etmodule = EvapotranspirationCO2(day, kiosk, parvalues)
        else:
            self.etmodule = Evapotranspiration(day, kiosk, parvalues)

    def __call__(self, day, drv):
        r = self.etmodule(day, drv)
        return r


class Evapotranspiration(SimulationObject):
    """Calculation of potential evaporation (water and soil) rates and actual
    crop transpiration rate.
        
    *Simulation parameters*:
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    CFET     Correction factor for potential transpiration   SCr       -
             rate.
    DEPNR    Dependency number for crop sensitivity to       SCr       -
             soil moisture stress.
    KDIFTB   Extinction coefficient for diffuse visible      TCr       -
             as function of DVS.
    IOX      Switch oxygen stress on (1) or off (0)          SCr       -
    IAIRDU   Switch airducts on (1) or off (0)               SCr       -
    CRAIRC   Critical air content for root aeration          SSo       -
    SM0      Soil porosity                                   SSo       -
    SMW      Volumetric soil moisture content at wilting     SSo       -
             point
    SMCFC    Volumetric soil moisture content at field       SSo       -
             capacity
    SM0      Soil porosity                                   SSo       -
    =======  ============================================= =======  ============
    

    *State variables*
    
    Note that these state variables are only assigned after finalize() has been
    run.

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    IDWST     Nr of days with water stress.                      N    -
    IDOST     Nr of days with oxygen stress.                     N    -
    =======  ================================================= ==== ============
    
    
    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    EVWMX    Maximum evaporation rate from an open water        Y    |cm day-1|
             surface.
    EVSMX    Maximum evaporation rate from a wet soil surface.  Y    |cm day-1|
    TRAMX    Maximum transpiration rate from the plant canopy   Y    |cm day-1|
    TRA      Actual transpiration rate from the plant canopy    Y    |cm day-1|
    IDOS     Indicates oxygen stress on this day (True|False)   N     -
    IDWS     Indicates water stress on this day (True|False)    N     -
    RFWS     Reduction factor for water stress                  N     -
    RFOS     Reduction factor for oxygen stress                 N     -
    RFTRA    Reduction factor for transpiration (wat & ox)      Y     -
    =======  ================================================= ==== ============
    
    *Signals send or handled*
    
    None
    
    *External dependencies:*
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    LAI      Leaf area index                     Leaf_dynamics       -
    SM       Volumetric soil moisture content    Waterbalance        -
    =======  =================================== =================  ============
    """

    # helper variable for Counting total days with water and oxygen
    # stress (IDWST, IDOST)
    _IDWST = Int(0)
    _IDOST = Int(0)

    class Parameters(ParamTemplate):
        CFET = Float(-99.)
        DEPNR = Float(-99.)
        KDIFTB = AfgenTrait()
        IAIRDU = Float(-99.)
        IOX = Float(-99.)
        CRAIRC = Float(-99.)
        SM0 = Float(-99.)
        SMW = Float(-99.)
        SMFCF = Float(-99.)

    class RateVariables(RatesTemplate):
        EVWMX = Float(-99.)
        EVSMX = Float(-99.)
        TRAMX = Float(-99.)
        TRA = Float(-99.)
        IDOS = Bool(False)
        IDWS = Bool(False)
        RFWS = Float(-99.)
        RFOS = Float(-99.)
        RFTRA = Float(-99.)

    class StateVariables(StatesTemplate):
        IDOST = Int(-99)
        IDWST = Int(-99)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["EVWMX", "EVSMX", "TRAMX", "TRA", "RFTRA"])
        self.states = self.StateVariables(kiosk, IDOST=-999, IDWST=-999)

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk
        
        KGLOB = 0.75 * p.KDIFTB(k.DVS)
  
        # crop specific correction on potential transpiration rate
        ET0_CROP = max(0., p.CFET * drv.ET0)
  
        # maximum evaporation and transpiration rates
        EKL = exp(-KGLOB * k.LAI)
        r.EVWMX = drv.E0 * EKL
        r.EVSMX = max(0., drv.ES0 * EKL)
        r.TRAMX = ET0_CROP * (1.-EKL)
                
        # Critical soil moisture
        SWDEP = SWEAF(ET0_CROP, p.DEPNR)
        SMCR = (1.-SWDEP)*(p.SMFCF-p.SMW) + p.SMW

        # Reduction factor for transpiration in case of water shortage (RFWS)
        r.RFWS = limit(0., 1., (k.SM - p.SMW)/(SMCR - p.SMW))

        # reduction in transpiration in case of oxygen shortage (RFOS)
        # for non-rice crops, and possibly deficient land drainage
        r.RFOS = 1.
        if p.IAIRDU == 0 and p.IOX == 1:
            RFOSMX = limit(0., 1., (p.SM0 - k.SM)/p.CRAIRC)
            # maximum reduction reached after 4 days
            r.RFOS = RFOSMX + (1. - min(k.DSOS, 4)/4.)*(1.-RFOSMX)

        # Transpiration rate multiplied with reduction factors for oxygen and water
        r.RFTRA = r.RFOS * r.RFWS
        r.TRA = r.TRAMX * r.RFTRA

        # Counting stress days
        if r.RFWS < 1.:
            r.IDWS = True
            self._IDWST += 1
        if r.RFOS < 1.:
            r.IDOS = True
            self._IDOST += 1

        return r.TRA, r.TRAMX
        
    @prepare_states
    def finalize(self, day):

        self.states.IDWST = self._IDWST
        self.states.IDOST = self._IDOST
        
        SimulationObject.finalize(self, day)


class EvapotranspirationCO2(SimulationObject):
    """Calculation of evaporation (water and soil) and transpiration rates
    taking into account the CO2 effect on crop transpiration.

    *Simulation parameters* (To be provided in cropdata dictionary):

    ======== ============================================= =======  ============
     Name     Description                                   Type     Unit
    ======== ============================================= =======  ============
    CFET     Correction factor for potential transpiration   S       -
             rate.
    DEPNR    Dependency number for crop sensitivity to       S       -
             soil moisture stress.
    KDIFTB   Extinction coefficient for diffuse visible      T       -
             as function of DVS.
    IOX      Switch oxygen stress on (1) or off (0)          S       -
    IAIRDU   Switch airducts on (1) or off (0)               S       -
    CRAIRC   Critical air content for root aeration          S       -
    SM0      Soil porosity                                   S       -
    SMW      Volumetric soil moisture content at wilting     S       -
             point
    SMCFC    Volumetric soil moisture content at field       S       -
             capacity
    SM0      Soil porosity                                   S       -
    CO2      Atmospheric CO2 concentration                   S       ppm
    CO2TRATB Reduction factor for TRAMX as function of
             atmospheric CO2 concentration                   T       -
    ======== ============================================= =======  ============


    *State variables*

    Note that these state variables are only assigned after finalize() has been
    run.

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    IDWST     Nr of days with water stress.                      N    -
    IDOST     Nr of days with oxygen stress.                     N    -
    =======  ================================================= ==== ============


    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    EVWMX    Maximum evaporation rate from an open water        Y    |cm day-1|
             surface.
    EVSMX    Maximum evaporation rate from a wet soil surface.  Y    |cm day-1|
    TRAMX    Maximum transpiration rate from the plant canopy   Y    |cm day-1|
    TRA      Actual transpiration rate from the plant canopy    Y    |cm day-1|
    IDOS     Indicates water stress on this day (True|False)    N    -
    IDWS     Indicates oxygen stress on this day (True|False)   N    -
    RFWS     Reducation factor for water stress                 Y     -
    RFOS     Reducation factor for oxygen stress                Y     -
    RFTRA    Reduction factor for transpiration (wat & ox)      Y     -
    =======  ================================================= ==== ============

    *Signals send or handled*

    None

    *External dependencies:*

    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    LAI      Leaf area index                     Leaf_dynamics       -
    SM       Volumetric soil moisture content    Waterbalance        -
    =======  =================================== =================  ============
    """

    # helper variable for counting total days with water and oxygen
    # stress (IDWST, IDOST)
    _IDWST = Int(0)
    _IDOST = Int(0)

    class Parameters(ParamTemplate):
        CFET    = Float(-99.)
        DEPNR   = Float(-99.)
        KDIFTB  = AfgenTrait()
        IAIRDU  = Float(-99.)
        IOX     = Float(-99.)
        CRAIRC  = Float(-99.)
        SM0     = Float(-99.)
        SMW     = Float(-99.)
        SMFCF   = Float(-99.)
        CO2     = Float(-99.)
        CO2TRATB = AfgenTrait()

    class RateVariables(RatesTemplate):
        EVWMX = Float(-99.)
        EVSMX = Float(-99.)
        TRAMX = Float(-99.)
        TRA   = Float(-99.)
        TRALY = Instance(array.array)
        IDOS  = Bool(False)
        IDWS  = Bool(False)
        RFWS = Float(-99.)
        RFOS = Float(-99.)
        RFTRA = Float(-99.)

    class StateVariables(StatesTemplate):
        IDOST  = Int(-99)
        IDWST  = Int(-99)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :param soildata: dictionary with WOFOST soildata key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["EVWMX","EVSMX", "TRAMX","TRA","TRALY", "RFTRA"])
        self.states = self.StateVariables(kiosk, IDOST=-999, IDWST=-999)

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk

        # reduction factor for CO2 on TRAMX
        RF_TRAMX_CO2 = p.CO2TRATB(p.CO2)

        # crop specific correction on potential transpiration rate
        ET0_CROP = max(0., p.CFET * drv.ET0)

        # maximum evaporation and transpiration rates
        KGLOB = 0.75*p.KDIFTB(k.DVS)
        EKL = exp(-KGLOB * k.LAI)
        r.EVWMX = drv.E0 * EKL
        r.EVSMX = max(0., drv.ES0 * EKL)
        r.TRAMX = ET0_CROP * (1.-EKL) * RF_TRAMX_CO2

        # Critical soil moisture
        SWDEP = SWEAF(ET0_CROP, p.DEPNR)

        SMCR = (1.-SWDEP)*(p.SMFCF-p.SMW) + p.SMW

        # Reduction factor for transpiration in case of water shortage (RFWS)
        r.RFWS = limit(0., 1., (k.SM-p.SMW)/(SMCR-p.SMW))

        # reduction in transpiration in case of oxygen shortage (RFOS)
        # for non-rice crops, and possibly deficient land drainage
        r.RFOS = 1.
        if p.IAIRDU == 0 and p.IOX == 1:
            RFOSMX = limit(0., 1., (p.SM0 - k.SM)/p.CRAIRC)
            # maximum reduction reached after 4 days
            r.RFOS = RFOSMX + (1. - min(k.DSOS, 4)/4.)*(1.-RFOSMX)

        # Transpiration rate multiplied with reduction factors for oxygen and water
        r.RFTRA = r.RFOS * r.RFWS
        r.TRA = r.TRAMX * r.RFTRA

        # Counting stress days
        if r.RFWS < 1.:
            r.IDWS = True
            self._IDWST += 1
        if r.RFOS < 1.:
            r.IDOS = True
            self._IDOST += 1

        return r.TRA, r.TRAMX

    @prepare_states
    def finalize(self, day):

        self.states.IDWST = self._IDWST
        self.states.IDOST = self._IDOST

        SimulationObject.finalize(self, day)


class EvapotranspirationCO2Layered(SimulationObject):
    """Calculation of evaporation (water and soil) and transpiration rates
    taking into account the CO2 effect on crop transpiration for a layered soil

    *Simulation parameters* (To be provided in cropdata dictionary):

    ======== ============================================= =======  ============
     Name     Description                                   Type     Unit
    ======== ============================================= =======  ============
    CFET     Correction factor for potential transpiration   S       -
             rate.
    DEPNR    Dependency number for crop sensitivity to       S       -
             soil moisture stress.
    KDIFTB   Extinction coefficient for diffuse visible      T       -
             as function of DVS.
    IOX      Switch oxygen stress on (1) or off (0)          S       -
    IAIRDU   Switch airducts on (1) or off (0)               S       -
    CRAIRC   Critical air content for root aeration          S       -
    SM0      Soil porosity                                   S       -
    SMW      Volumetric soil moisture content at wilting     S       -
             point
    SMCFC    Volumetric soil moisture content at field       S       -
             capacity
    SM0      Soil porosity                                   S       -
    CO2      Atmospheric CO2 concentration                   S       ppm
    CO2TRATB Reduction factor for TRAMX as function of
             atmospheric CO2 concentration                   T       -
    ======== ============================================= =======  ============


    *State variables*

    Note that these state variables are only assigned after finalize() has been
    run.

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    IDWST     Nr of days with water stress.                      N    -
    IDOST     Nr of days with oxygen stress.                     N    -
    =======  ================================================= ==== ============


    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    EVWMX    Maximum evaporation rate from an open water        Y    |cm day-1|
             surface.
    EVSMX    Maximum evaporation rate from a wet soil surface.  Y    |cm day-1|
    TRAMX    Maximum transpiration rate from the plant canopy   Y    |cm day-1|
    TRA      Actual transpiration rate from the plant canopy    Y    |cm day-1|
    IDOS     Indicates water stress on this day (True|False)    N    -
    IDWS     Indicates oxygen stress on this day (True|False)   N    -
    RFWS     Reducation factor for water stress                 Y     -
    RFOS     Reducation factor for oxygen stress                Y     -
    RFTRA    Reduction factor for transpiration (wat & ox)      Y     -
    =======  ================================================= ==== ============

    *Signals send or handled*

    None

    *External dependencies:*

    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    LAI      Leaf area index                     Leaf_dynamics       -
    SM       Volumetric soil moisture content    Waterbalance        -
    =======  =================================== =================  ============
    """

    # helper variable for counting total days with water and oxygen
    # stress (IDWST, IDOST)
    _IDWST = Int(0)
    _IDOST = Int(0)
    _DSOS = Int(0)

    soil_profile = None

    class Parameters(ParamTemplate):
        CFET    = Float(-99.)
        DEPNR   = Float(-99.)
        KDIFTB  = AfgenTrait()
        IAIRDU  = Float(-99.)
        IOX     = Float(-99.)
        CO2     = Float(-99.)
        CO2TRATB = AfgenTrait()

    class RateVariables(RatesTemplate):
        EVWMX = Float(-99.)
        EVSMX = Float(-99.)
        TRAMX = Float(-99.)
        TRA   = Float(-99.)
        TRALY = Instance(np.ndarray)
        IDOS  = Bool(False)
        IDWS  = Bool(False)
        RFWS = Instance(np.ndarray)
        RFOS = Instance(np.ndarray)
        RFTRALY = Instance(np.ndarray)
        RFTRA = Float(-99.)

    class StateVariables(StatesTemplate):
        IDOST  = Int(-99)
        IDWST  = Int(-99)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :param soildata: dictionary with WOFOST soildata key/value pairs
        """

        self.soil_profile = parvalues["soil_profile"]
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["EVWMX","EVSMX", "TRAMX","TRA","TRALY", "RFTRA"])
        self.states = self.StateVariables(kiosk, IDOST=-999, IDWST=-999)

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk

        # reduction factor for CO2 on TRAMX
        RF_TRAMX_CO2 = p.CO2TRATB(p.CO2)

        # crop specific correction on potential transpiration rate
        ET0_CROP = max(0., p.CFET * drv.ET0)

        # maximum evaporation and transpiration rates
        KGLOB = 0.75*p.KDIFTB(k.DVS)
        EKL = exp(-KGLOB * k.LAI)
        r.EVWMX = drv.E0 * EKL
        r.EVSMX = max(0., drv.ES0 * EKL)
        r.TRAMX = ET0_CROP * (1.-EKL) * RF_TRAMX_CO2

        # Critical soil moisture
        SWDEP = SWEAF(ET0_CROP, p.DEPNR)
        depth = 0.0

        RFWS = np.zeros(len(self.soil_profile), dtype=np.float64)
        RFOS = np.zeros_like(RFWS)
        TRALY = np.zeros_like(RFWS)

        layercnt = range(len(self.soil_profile))
        for i, SM, layer in zip(layercnt, k.SM, self.soil_profile):
            SMCR = (1.-SWDEP)*(layer.SMFCF - layer.SMW) + layer.SMW

            # Reduction factor for transpiration in case of water shortage (RFWS)
            RFWS[i] = limit(0., 1., (SM - layer.SMW)/(SMCR - layer.SMW))

            # reduction in transpiration in case of oxygen shortage (RFOS)
            # for non-rice crops, and possibly deficient land drainage
            RFOS[i] = 1.
            if p.IAIRDU == 0 and p.IOX == 1:
                SMAIR = layer.SM0 - layer.CRAIRC
                if(SM >= SMAIR):
                    self._DSOS = min((self._DSOS + 1), 4)
                else:
                    self._DSOS = 0                
                RFOSMX = limit(0., 1., (layer.SM0 - SM)/layer.CRAIRC)
                RFOS[i] = RFOSMX + (1. - min( self._DSOS, 4)/4.)*(1.-RFOSMX)
            root_fraction = max(0.0, (min(k.RD, depth + layer.Thickness) - depth)) / k.RD
            RFTRA_layer = RFOS[i] * RFWS[i]
            TRALY[i] = r.TRAMX * RFTRA_layer * root_fraction

            depth += layer.Thickness

        r.TRA = TRALY.sum()
        r.TRALY = TRALY
        r.RFTRA = r.TRA/r.TRAMX if r.TRAMX > 0. else 1.
        r.RFOS = RFOS
        r.RFWS = RFWS

        # Counting stress days
        if any(r.RFWS < 1.):
            r.IDWS = True
            self._IDWST += 1
        if any(r.RFOS < 1.):
            r.IDOS = True
            self._IDOST += 1

    @prepare_states
    def finalize(self, day):

        self.states.IDWST = self._IDWST
        self.states.IDOST = self._IDOST

        SimulationObject.finalize(self, day)


class Simple_Evapotranspiration(SimulationObject):
    """Calculation of evaporation (water and soil) and transpiration rates.
    
    simplified compared to the WOFOST model such that the parameters are not
    needed to calculate the ET. Parameters such as KDIF, CFET, DEPNR have been
    hardcoded and taken from a typical cereal crop. Also the oxygen stress has
    been switched off.
    
    *Simulation parameters* (To be provided in soildata dictionary):
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    SMW      Volumetric soil moisture content at wilting     S       -
             point
    SMCFC    Volumetric soil moisture content at field       S       -
             capacity
    SM0      Soil porosity                                   S       -
    =======  ============================================= =======  ============
    

    *State variables*
    
    None
    
    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    EVWMX    Maximum evaporation rate from an open water        Y    |cm day-1|
             surface.
    EVSMX    Maximum evaporation rate from a wet soil surface.  Y    |cm day-1|
    TRAMX    Maximum transpiration rate from the plant canopy   Y    |cm day-1|
    TRA      Actual transpiration rate from the plant canopy    Y    |cm day-1|
    =======  ================================================= ==== ============
    
    *Signals send or handled*
    
    None
    
    *External dependencies:*
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    LAI      Leaf area index                     Leaf_dynamics       -
    SM       Volumetric soil moisture content    Waterbalance        -
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):
        SM0    = Float(-99.)
        SMW    = Float(-99.)
        SMFCF  = Float(-99.)
        CFET = Float(-99.)
        DEPNR = Float(-99.)

    class RateVariables(RatesTemplate):
        EVWMX = Float(-99.)
        EVSMX = Float(-99.)
        TRAMX = Float(-99.)
        TRA   = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param soildata: dictionary with WOFOST soildata key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["EVWMX","EVSMX",
                                                        "TRAMX","TRA"])

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        
        LAI = self.kiosk["LAI"]
        SM  = self.kiosk["SM"]
        
        # value for KDIF taken for maize
        KDIF = 0.6
        KGLOB = 0.75*KDIF
  
        # crop specific correction on potential transpiration rate
        ET0 = p.CFET * drv.ET0
  
        # maximum evaporation and transpiration rates
        EKL = exp(-KGLOB * LAI)
        r.EVWMX = drv.E0 * EKL
        r.EVSMX = max(0., drv.ES0 * EKL)
        r.TRAMX = max(0.000001, ET0 * (1.-EKL))
        
        # Critical soil moisture
        SWDEP = SWEAF(ET0, p.DEPNR)
        
        SMCR = (1.-SWDEP)*(p.SMFCF-p.SMW) + p.SMW

        # Reduction factor for transpiration in case of water shortage (RFWS)
        RFWS = limit(0., 1., (SM-p.SMW)/(SMCR-p.SMW))
        
        # Reduction factor for oxygen stress set to 1.
        RFOS = 1.0

        # Transpiration rate multiplied with reduction factors for oxygen and
        # water
        r.TRA = r.TRAMX * RFOS * RFWS

        return (r.TRA, r.TRAMX)