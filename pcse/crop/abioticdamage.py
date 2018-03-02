# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""
Components for modelling of abiotic damage to crops. 

The following components are available:
* Frost damage:
  - FROSTOL: models LT50 to estimate leaf and crop death 
  - CERES_WinterKill: models a hardening index to estimate leaf and crop death
"""

#!/usr/bin/env python
import os
from math import exp

from ..traitlets import Float, Int, Instance, Enum, Bool
from ..decorators import prepare_rates, prepare_states

from ..util import limit, merge_dict
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject, VariableKiosk
from .. import signals
from .. import exceptions as exc

class CrownTemperature(SimulationObject):
    """Implementation of a simple algorithm for estimating the crown temperature
    (2cm under the soil surface) under snow.
    
    Is is based on a simple empirical equation which estimates the daily
    minimum, maximum and mean crown 
    temperature as a function of daily min or max temperature and the relative
    snow depth (RSD):
    
    :math:`RSD = min(15, SD)/15`
    
    and
    
    :math:`T^{crown}_{min} = T_{min} * (A + B(1 - RSD)^{2})`   

    and

    :math:`T^{crown}_{max} = T_{max} * (A + B(1 - RSD)^{2})`
    
    and
    
    :math:`T^{crown}_{avg} = (T^{crown}_{max} + T^{crown}_{min})/2`

    At zero snow depth crown temperature is estimated close the the air
    temperature. Increasing snow depth acts as a buffer damping the effect of
    low air temperature on the crown temperature. The maximum value of the
    snow depth is limited on 15cm. Typical values for A and B are 0.2 and
    0.5

    Note that the crown temperature is only estimated if drv.TMIN<0, otherwise
    the TMIN, TMAX and daily average temperature (TEMP) are returned.

    :param day: day when model is initialized
    :param kiosk: VariableKiosk of this instance
    :param parvalues: `ParameterProvider` object providing parameters as
            key/value pairs
    :returns: a tuple containing minimum, maximum and daily average crown
              temperature. 

    *Simulation parameters*
    
    ========= ============================================== =======  ==========
     Name     Description                                    Type     Unit
    ========= ============================================== =======  ==========
    ISNOWSRC  Use prescribed snow depth from driving          SSi      -
              variables (0) or modelled snow depth through
              the kiosk (1)
    CROWNTMPA A parameter in equation for crown temperature   SSi      -
    CROWNTMPB B parameter in equation for crown temperature   SSi      -
    ========= ============================================== =======  ==========

    *Rate variables*

    ========== =============================================== =======  ==========
     Name      Description                                      Pbl     Unit
    ========== =============================================== =======  ==========
    TEMP_CROWN  Daily average crown temperature                   N       |C|
    TMIN_CROWN  Daily minimum crown temperature                   N       |C|
    TMAX_CROWN  Daily maximum crown temperature                   N       |C|
    ========== =============================================== =======  ==========

    Note that the calculated crown temperatures are not real rate variables as
    they do not pertain to rate of change. In fact they are a `derived driving
    variable`. Nevertheless for calculating the frost damage they should
    become available during the rate calculation step and by treating them
    as rate variables, they can be found by a `get_variable()` call and thus
    be defined in the list of OUTPUT_VARS in the configuration file

    *External dependencies:*

    ============ =============================== ========================== =====
     Name        Description                         Provided by             Unit
    ============ =============================== ========================== =====
    SNOWDEPTH    Depth of snow cover.             Prescibed by driving       |cm|
                                                  variables or simulated
                                                  by snow cover module and
                                                  taken from kiosk
    ============ =============================== ========================== =====
    """
    # This setting is only used when running the unit tests for FROSTOL.
    # For unit testing, FROSTOL should not rely on the CrownTemperature model,
    # but instead use the prescribed crown temperature directly.
    _testing_ = Bool(False)

    class Parameters(ParamTemplate):
        CROWNTMPA = Float()
        CROWNTMPB = Float()
        ISNOWSRC  = Float()

    class RateVariables(RatesTemplate):
        TEMP_CROWN = Float()
        TMIN_CROWN = Float()
        TMAX_CROWN = Float()

    def initialize(self, day, kiosk, parvalues, testing=False):
        self.kiosk = kiosk
        self._testing_ = testing
        if not self._testing_:
            self.params = self.Parameters(parvalues)
            self.rates = self.RateVariables(self.kiosk)

    @prepare_rates
    def __call__(self, day, drv):

        # If unit testing then directly return the prescribed crown temperature
        if self._testing_:
            return (0., 10., drv.TEMP_CROWN)

        p = self.params
        r = self.rates

        # Take snow depth from driving variables or kiosk depending on
        # ISNOWSRC and limit the snow depth on 15 cm
        if p.ISNOWSRC == 0:
            SD = drv.SNOWDEPTH
        else:
            SD = self.kiosk["SNOWDEPTH"]
        RSD = limit(0., 15., SD)/15.

        if drv.TMIN < 0:
            r.TMIN_CROWN = drv.TMIN*(p.CROWNTMPA + p.CROWNTMPB*(1. - RSD)**2)
            r.TMAX_CROWN = drv.TMAX*(p.CROWNTMPA + p.CROWNTMPB*(1. - RSD)**2)
            r.TEMP_CROWN = (r.TMIN_CROWN + r.TMAX_CROWN)/2.
        else:
            r.TMIN_CROWN = drv.TMIN
            r.TMAX_CROWN = drv.TMAX
            r.TEMP_CROWN = drv.TEMP

        return (r.TMIN_CROWN, r.TMAX_CROWN, r.TEMP_CROWN)


class FROSTOL(SimulationObject):
    """ Implementation of the FROSTOL model for frost damage in winter-wheat.

    :param day: start date of the simulation
    :param kiosk: variable kiosk of this PCSE instance
    :param parvalues: `ParameterProvider` object providing parameters as
            key/value pairs

    *Simulation parameters*
    
    ============== ============================================= =======  ============
     Name          Description                                   Type     Unit
    ============== ============================================= =======  ============
    IDSL           Switch for phenological development options    SCr      -
                   temperature only (IDSL=0), including           
                   daylength (IDSL=1) and including               
                   vernalization (IDSL>=2). FROSTOL requires
                   IDSL>=2
    LT50C          Critical LT50 defined as the lowest LT50       SCr     |C|
                   value that the wheat cultivar can obtain
    FROSTOL_H      Hardening coefficient                          SCr     |C-1day-1| 
    FROSTOL_D      Dehardening coefficient                        SCr     |C-3day-1|
    FROSTOL_S      Low temperature stress coefficient             SCr     |C-1day-1|
    FROSTOL_R      Respiration stress coefficient                 SCr     |day-1|
    FROSTOL_SDBASE Minimum snow depth for respiration stress      SCr      cm
    FROSTOL_SDMAX  Snow depth with maximum respiration stress.    SCr      cm
                   Larger snow depth does not increase stress
                   anymore.
    FROSTOL_KILLCF Steepness coefficient for logistic kill        SCr     -
                   function.
    ISNOWSRC       Use prescribed snow depth from driving         SSi     -
                   variables (0) or modelled snow depth through
                   the kiosk (1)
    ============== ============================================= =======  ============

    *State variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
     LT50T    Current LT50 value                                N    |C|
     LT50I    Initial LT50 value of unhardened crop             N    |C|
     IDFST    Total number of days with frost stress            N    -
    =======  ================================================= ==== ============


    *Rate variables*

    ========== ================================================= ==== ============
     Name       Description                                      Pbl      Unit
    ========== ================================================= ==== ============
    RH         Rate of hardening                                  N    |C day-1|
    RDH_TEMP   Rate of dehardening due to temperature             N    |C day-1|
    RDH_RESP   Rate of dehardening due to respiration stress      N    |C day-1|
    RDH_TSTR   Rate of dehardening due to temperature stress      N    |C day-1|
    IDFS       Frost stress, yes (1) or no (0). Frost stress is   N    -
               defined as: RF_FROST > 0
    RF_FROST   Reduction factor on leave biomass as a function    Y    -
               of min. crown temperature and LT50T: ranges
               from 0 (no damage) to 1 (complete kill).
    RF_FROST_T Total frost kill through the growing season        N    -
               is computed as the multiplication of the daily
               frost kill events, 0 means no damage, 1 means
               total frost kill.
    ========== ================================================= ==== ============

    
    *External dependencies:*
    
    ============ =============================== ========================== =====
     Name        Description                         Provided by             Unit
    ============ =============================== ========================== =====
    TEMP_CROWN   Daily average crown temperature  CrownTemperature           |C|
                 derived from calling the
                 crown_temperature module.
    TMIN_CROWN   Daily minimum crown temperature  CrownTemperature           |C|
                 derived from calling the
                 crown_temperature module.
    ISVERNALISED Boolean reflecting the
                 vernalisation state of the       Vernalisation i.c.m. with   -
                 crop.                            DVS_Phenology module
    ============ =============================== ========================== =====

    Reference: Anne Kari Bergjord, Helge Bonesmo, Arne Oddvar Skjelvag, 2008.
               Modelling the course of frost tolerance in winter wheat: I. Model
               development, European Journal of Agronomy, Volume 28,
               Issue 3, April 2008, Pages 321-330.
    
    http://dx.doi.org/10.1016/j.eja.2007.10.002
    """
    crown_temperature = Instance(SimulationObject)

    # Helper variable for remaining crop fraction as a result of frost kill
    _CROP_FRACTION_REMAINING = Float(1.0)

    class Parameters(ParamTemplate):
        IDSL      = Float(-99.)
        LT50C     = Float(-99.)
        FROSTOL_H = Float(-99.)
        FROSTOL_D = Float(-99.)
        FROSTOL_S = Float(-99.)
        FROSTOL_R = Float(-99.)
        FROSTOL_SDBASE = Float(-99.)
        FROSTOL_SDMAX  = Float(-99.)
        FROSTOL_KILLCF = Float(-99)
        ISNOWSRC = Float(-99)

    class RateVariables(RatesTemplate):
        RH       = Float(-99.)
        RDH_TEMP = Float(-99.)
        RDH_RESP = Float(-99.)
        RDH_TSTR = Float(-99.)
        IDFS     = Int(-99)
        RF_FROST = Float(-99.)

    class StateVariables(StatesTemplate):
        LT50T = Float(-99.)
        LT50I = Float(-99.)
        IDFST = Int(-99)
        RF_FROST_T = Float(-99)

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues, testing=False):

        # Initialize module for crown temperature
        self.crown_temperature = CrownTemperature(day, kiosk, parvalues,
                                                  testing)

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="RF_FROST")
        self.kiosk = kiosk

        # Define initial states
        LT50I = -0.6 + 0.142 * self.params.LT50C
        self.states = self.StateVariables(kiosk, LT50T=LT50I, LT50I=LT50I,
                                          IDFST=0, RF_FROST_T=0.)

        # Check on vernalization
        if self.params.IDSL < 2:
            msg = ("FROSTOL needs vernalization to be enabled in the " +
                   "phenology module (IDSL=2).")
            self.logger.error(msg)
            raise exc.ParameterError(msg)

    #---------------------------------------------------------------------------
    @prepare_rates
    def calc_rates(self, day, drv):

        r = self.rates
        p = self.params
        s = self.states

        # vernalisation state
        isVernalized = self.kiosk["ISVERNALISED"]

        TMIN_CROWN, TMAX_CROWN, TEMP_CROWN = self.crown_temperature(day, drv)

        # p.ISNOWSRC=0 derive snow depth from driving variables `drv`
        # else assume snow depth is a published state variable
        if p.ISNOWSRC == 0:
            snow_depth = drv.SNOWDEPTH
        else:
            snow_depth = self.kiosk["SNOWDEPTH"]

        # Hardening
        if (not isVernalized) and (TEMP_CROWN < 10.):
            xTC = limit(0., 10., TEMP_CROWN)
            r.RH = p.FROSTOL_H * (10. - xTC)*(s.LT50T - p.LT50C)
        else:
            r.RH = 0.

        # Dehardening
        TCcrit = (10. if (not isVernalized) else -4.)
        if TEMP_CROWN > TCcrit:
            r.RDH_TEMP = p.FROSTOL_D * (s.LT50I - s.LT50T) * \
                         (TEMP_CROWN + 4)**3
        else:
            r.RDH_TEMP = 0.

        # Stress due to respiration under snow coverage
        xTC = (TEMP_CROWN if TEMP_CROWN > -2.5 else -2.5)
        Resp = (exp(0.84 + 0.051*xTC)-2.)/1.85

        Fsnow = (snow_depth - p.FROSTOL_SDBASE)/(p.FROSTOL_SDMAX - p.FROSTOL_SDBASE)
        Fsnow = limit(0., 1., Fsnow)
        r.RDH_RESP = p.FROSTOL_R * Resp * Fsnow

        # Stress due to low temperatures
        r.RDH_TSTR = (s.LT50T - TEMP_CROWN) * \
                      1./exp(-p.FROSTOL_S * (s.LT50T - TEMP_CROWN) - 3.74)

        # kill factor using logistic function. Because the logistic function
        # stretches from -inf to inf, some limits must be applied. In this
        # case we assume that killfactor < 0.05 means no kill and
        # killfactor > 0.95 means complete kill.
        if TMIN_CROWN < 0.:
            killfactor = 1/(1 + exp((TMIN_CROWN - s.LT50T)/p.FROSTOL_KILLCF))
            if killfactor < 0.05:
                killfactor = 0.
            elif killfactor > 0.95:
                killfactor = 1.
        else:
            killfactor = 0.

        # Frost stress occurring yes/no
        r.IDFS = 1 if (killfactor > 0.) else 0

        # Reduction factor on leave biomass
        r.RF_FROST = killfactor

        # Fraction of the remaining standing crop
        self._CROP_FRACTION_REMAINING *= (1. - killfactor)

    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day, delt=1.0):
        states = self.states
        rates  = self.rates
        params = self.params

        # Change hardening state
        LT50T = states.LT50T
        LT50T -= rates.RH
        LT50T += (rates.RDH_TEMP + rates.RDH_RESP + rates.RDH_TSTR)
        states.LT50T = limit(params.LT50C, states.LT50I, LT50T)

        # Count number of days with frost stress
        states.IDFST += rates.IDFS

        # Total cumulative frost kill computed as 1 minus the fraction
        # of remaining living crop
        states.RF_FROST_T = 1. - self._CROP_FRACTION_REMAINING


class CERES_WinterKill(SimulationObject):
    """Implementation of the winter-kill module in the CERES-wheat model (CWWK).

    :param day: start date of the simulation
    :param kiosk: variable kiosk of this PCSE instance
    :param parvalues: `ParameterProvider` object providing parameters as
            key/value pairs

    *Simulation parameters*

    ============== ============================================= =======  ============
     Name          Description                                   Type     Unit
    ============== ============================================= =======  ============
    CWWK_HC_S1      CERES hardening coefficient for stage 1        SCr      TBD
    CWWK_HC_S2      CERES hardening coefficient for stage 1        SCr      TBD
    CWWK_DHC        CERES dehardening coefficient                  Scr      TBD
    CWWK_KILLTEMP   CERES Killing temperature per HI               Scr      |C|
    ============== ============================================= =======  ============

    *State variables*

    =========== ================================================= ==== ======
     Name          Description                                     Pbl  Unit
    =========== ================================================= ==== ======
     HARDINDEX    Hardening index                                   N     -
     HIKILLTEMP   Killing temperature as function of HI             N    |C|
    =========== ================================================= ==== ======


    *Rate variables*

    ============ ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    ============ ================================================= ==== ============
    RH           Rate of hardening                                  N    |day-1|
    RDH          Rate of dehardening                                N    |day-1|
    HIKILLFACTOR Fraction of biomass killed by low temp             N    -
    ============ ================================================= ==== ============

    Reference:
    Savdie, I., R. Whitewood, et al. (1991). Potential for winter wheat 
    production in western Canada: A CERES model winterkill risk 
    assessment. Canadian Journal of Plant Science 71: 21-30.
    """

    class Parameters(ParamTemplate):
        CWWK_HC_S1  = Float(-99.) # Hardening coefficient stage 1
        CWWK_HC_S2  = Float(-99.) # Hardening coefficient stage 2
        CWWK_DHC = Float(-99.)    # De-hardening coefficient
        CWWK_KILLTEMP = Float(-99.) # Initial Killing temperature

    class StateVariables(StatesTemplate):
        HARDINDEX  = Float(-99.) # Hardening Index
        HIKILLTEMP = Float(-99.) # Kill temperature given Hardening Index

    class RateVariables(RatesTemplate):
        RH = Float(-99.)
        RDH = Float(-99.)
        HIKILLFACTOR = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="HIKILLFACTOR")
        self.kiosk = kiosk

        # Define initial states
        self.states = self.StateVariables(kiosk, HARDINDEX=0.,
                                          HIKILLTEMP=self.params.CWWK_KILLTEMP)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        params = self.params
        states = self.states

        # derive snow depth from kiosk
        snow_depth = self.kiosk["SNOWDEPTH"]

        if states.HARDINDEX >= 1.: # HI between 1 and 2.
            if drv.TEMP_CROWN < 0.:
                # 12 days of hardening are enough to reach stage 2
                # default value 0.083333 = 1/12
                rates.RH = params.CWWK_HC_S2
            else:
                rates.RH = 0.
        else:  # HI between 0 and 1
            if (drv.TEMP_CROWN > -1.) and (drv.TEMP_CROWN < 8.):
                # At 3.5 degree HI increase 0.1 (max) and with 0.06 (min) 
                # at -1 and 8 degree. Default vaue for CERESWK_HC_S1=0.1
                rates.RH = params.CWWK_HC_S1 - \
                                       ((3.5 - drv.TEMP_CROWN)**2/506.)
            else:
                rates.RH = 0.

        # Dehardening
        if drv.TMAX_CROWN > 10:
            #for each degree above 10, HI decreases with 0.02
            rates.RDH = (10 - drv.TMAX_CROWN) * params.CWWK_DHC
        else:
            rates.RDH = 0.

        # Calculate the killing factor based on the current kill temperature
        if drv.TMIN_CROWN < states.HIKILLTEMP:
            rates.KILLFACTOR = 1.
            # Send signal that crop is finished
            self._send_signal(signals.crop_finish, day=day, finish_type="frost kill")

        elif drv.TMIN_CROWN > params.CWWK_KILLTEMP:
            rates.KILLFACTOR = 0.

        else:
            KF = (0.02 * states.HARDINDEX - 0.1) * \
                  ((drv.TMINCROWN * 0.85) + (drv.TMAX_CROWN * 0.15) + \
                   10 + (0.25 * snow_depth))
            rates.KILLFACTOR = limit(0, 0.96, KF)

    @prepare_states
    def integrate(self, day, delt=1.0):
        states = self.states
        rates  = self.rates
        params = self.params

        states.HARDINDEX += (rates.RH + rates.RDH)
        states.HIKILLTEMP = (states.HARDINDEX + 1.) * params.CWWK_KILLTEMP

