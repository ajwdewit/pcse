# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Python implementations of the WOFOST waterbalance modules for simulation
of potential production (`WaterbalancePP`) and water-limited production
(`WaterbalanceFD`) under freely draining conditions.
"""
from math import sqrt

from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, List
from ..decorators import prepare_rates, prepare_states
from ..util import limit, Afgen, merge_dict
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals
from .. import exceptions as exc
from .snowmaus import SnowMAUS


class WaterbalancePP(SimulationObject):
    """Fake waterbalance for simulation under potential production.

    Keeps the soil moisture content at field capacity and only accumulates crop transpiration
    and soil evaporation rates through the course of the simulation
    """
    # Counter for Days-Dince-Last-Rain
    DSLR = Float(1)
    # rainfall rate of previous day
    RAINold = Float(0)

    class Parameters(ParamTemplate):
        SMFCF = Float(-99.)

    class StateVariables(StatesTemplate):
        SM = Float(-99.)
        WTRAT = Float(-99.)
        EVST = Float(-99.)

    class RateVariables(RatesTemplate):
        EVS = Float(-99.)
        WTRA = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """    
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: ParameterProvider object containing all parameters
    
        This waterbalance keeps the soil moisture always at field capacity. Therefore   
        `WaterbalancePP` has only one parameter (`SMFCF`: the field capacity of the
        soil) and one state variable (`SM`: the volumetric soil moisture).
        """
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="EVS")
        self.states = self.StateVariables(kiosk, SM=self.params.SMFCF,
                                          publish="SM", EVST=0, WTRAT=0)
    
    @prepare_rates
    def calc_rates(self, day, drv):

        r = self.rates
        # Transpiration and maximum soil and surface water evaporation rates
        # are calculated by the crop Evapotranspiration module.
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if "TRA" not in self.kiosk:
            r.WTRA = 0.
            EVSMX = drv.ES0
        else:
            r.WTRA = self.kiosk["TRA"]
            EVSMX = self.kiosk["EVSMX"]

        # Actual evaporation rates

        if self.RAINold >= 1:
            # If rainfall amount >= 1cm on previous day assume maximum soil
            # evaporation
            r.EVS = EVSMX
            self.DSLR = 1.
        else:
            # Else soil evaporation is a function days-since-last-rain (DSLR)
            self.DSLR += 1
            EVSMXT = EVSMX * (sqrt(self.DSLR) - sqrt(self.DSLR - 1))
            r.EVS = min(EVSMX, EVSMXT + self.RAINold)

        # Hold rainfall amount to keep track of soil surface wetness and reset self.DSLR if needed
        self.RAINold = drv.RAIN

    @prepare_states
    def integrate(self, day, delt=1.0):

        # Keep soil moisture on field capacity
        self.states.SM = self.params.SMFCF

        # Accumulated transpiration and soil evaporation amounts
        self.states.EVST += self.rates.EVS * delt
        self.states.WTRAT += self.rates.WTRA * delt
        

class WaterbalanceFD(SimulationObject):
    """Waterbalance for freely draining soils under water-limited production.

    The purpose of the soil water balance calculations is to estimate the
    daily value of the soil moisture content. The soil moisture content
    influences soil moisture uptake and crop transpiration.
    
    The dynamic calculations are carried out in two sections, one for the 
    calculation of rates of change per timestep (= 1 day) and one for the
    calculation of summation variables and state variables. The water balance
    is driven by rainfall, possibly buffered as surface storage, and
    evapotranspiration. The processes considered are infiltration, soil water
    retention, percolation (here conceived as downward water flow from rooted
    zone to second layer), and the loss of water beyond the maximum root zone. 

    The textural profile of the soil is conceived as homogeneous. Initially the
    soil profile consists of two layers, the actually rooted soil and the soil
    immediately below the rooted zone until the maximum rooting depth is reached
    by roots(soil and crop dependent). The extension of the root zone from the
    initial rooting depth to maximum rooting depth is described in Root_Dynamics
    class. From the moment that the maximum rooting depth is reached the soil
    profile may be described as a one layer system depending if the roots are
    able to penetrate the entire profile. If not a non-rooted part remains
    at the bottom of the profile.

    The class WaterbalanceFD is derived from WATFD.FOR in WOFOST7.1 with the
    exception that the depth of the soil is now completely determined by the
    maximum soil depth (RDMSOL) and not by the minimum of soil depth and crop
    maximum rooting depth (RDMCR).
    
    **Simulation parameters:**
    
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    SMFCF     Field capacity of the soil                       SSo     -
    SM0       Porosity of the soil                             SSo     -
    SMW       Wilting point of the soil                        SSo     -
    CRAIRC    Soil critical air content (waterlogging)         SSo     -
    SOPE      maximum percolation rate root zone               SSo    |cmday-1|
    KSUB      maximum percolation rate subsoil                 SSo    |cmday-1|
    RDMSOL    Soil rootable depth                              SSo     cm
    IFUNRN    Indicates whether non-infiltrating fraction of   SSi    -
              rain is a function of storm size (1)
              or not (0)                                      
    SSMAX     Maximum surface storage                          SSi     cm
    SSI       Initial surface storage                          SSi     cm
    WAV       Initial amount of water in total soil            SSi     cm
              profile
    NOTINF    Maximum fraction of rain not-infiltrating into   SSi     -
              the soil
    SMLIM     Initial maximum moisture content in initial      SSi     -
              rooting depth zone.
    ======== =============================================== =======  ==========

    **State variables:**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    SM        Volumetric moisture content in root zone          Y    -
    SS        Surface storage (layer of water on surface)       N    cm
    SSI       Initial urface storage                            N    cm
    W         Amount of water in root zone                      N    cm
    WI        Initial amount of water in the root zone          N    cm
    WLOW      Amount of water in the subsoil (between current   N    cm
              rooting depth and maximum rootable depth)
    WLOWI     Initial amount of water in the subsoil                 cm
    WWLOW     Total amount of water in the  soil profile        N    cm
              WWLOW = WLOW + W
    WTRAT     Total water lost as transpiration as calculated   N    cm
              by the water balance. This can be different 
              from the CTRAT variable which only counts
              transpiration for a crop cycle.
    EVST      Total evaporation from the soil surface           N    cm
    EVWT      Total evaporation from a water surface            N    cm
    TSR       Total surface runoff                              N    cm
    RAINT     Total amount of rainfall (eff + non-eff)          N    cm
    WDRT      Amount of water added to root zone by increase    N    cm
              of root growth
    TOTINF    Total amount of infiltration                      N    cm
    TOTIRR    Total amount of effective irrigation              N    cm
    PERCT     Total amount of water percolating from rooted     N    cm
              zone to subsoil
    LOSST     Total amount of water lost to deeper soil         N    cm
    DSOS      Days since oxygen stress, accumulates the number  Y     -
              of consecutive days of oxygen stress
    WBALRT    Checksum for root zone waterbalance. Will be      N    cm
              calculated within `finalize()`, abs(WBALRT) >
              0.0001 will raise a WaterBalanceError
    WBALTT    Checksum for total waterbalance. Will be          N    cm
              calculated within `finalize()`, abs(WBALTT) >
              0.0001 will raise a WaterBalanceError
    =======  ================================================= ==== ============

    **Rate variables:**

    =========== ================================================= ==== ============
     Name        Description                                      Pbl      Unit
    =========== ================================================= ==== ============
    EVS         Actual evaporation rate from soil                  N    |cmday-1|
    EVW         Actual evaporation rate from water surface         N    |cmday-1|
    WTRA        Actual transpiration rate from plant canopy,       N    |cmday-1|
                is directly derived from the variable "TRA" in
                the evapotranspiration module
    RAIN_INF    Infiltrating rainfall rate for current day         N    |cmday-1|
    RAIN_NOTINF Non-infiltrating rainfall rate for current day   N    |cmday-1|
    RIN         Infiltration rate for current day                  N    |cmday-1|
    RIRR        Effective irrigation rate for current day,         N    |cmday-1|
                computed as irrigation amount * efficiency.
    PERC        Percolation rate to non-rooted zone                N    |cmday-1|
    LOSS        Rate of water loss to deeper soil                  N    |cmday-1|
    DW          Change in amount of water in rooted zone as a      N    |cmday-1|
                result of infiltration, transpiration and
                evaporation.
    DWLOW       Change in amount of water in subsoil               N    |cmday-1|
    DTSR        Change in surface runoff                           N    |cmday-1|
    DSS         Change in surface storage                          N    |cmday-1|
    =========== ================================================= ==== ============
    
    
    **External dependencies:**
    
    ============ ============================== ====================== =========
     Name        Description                         Provided by         Unit
    ============ ============================== ====================== =========
     TRA          Crop transpiration rate       Evapotranspiration     |cmday-1|
     EVSMX        Maximum evaporation rate      Evapotranspiration     |cmday-1|
                  from a soil surface below
                  the crop canopy
     EVWMX        Maximum evaporation rate       Evapotranspiration    |cmday-1|
                  from a water surface below
                  the crop canopy
     RD           Rooting depth                  Root_dynamics          cm
    ============ ============================== ====================== =========

    **Exceptions raised:**
    
    A WaterbalanceError is raised when the waterbalance is not closing at the
    end of the simulation cycle (e.g water has "leaked" away).
    """
    # previous and maximum rooting depth value
    RDold = Float(-99.)
    RDM = Float(-99.)
    # Counter for Days-Dince-Last-Rain 
    DSLR = Float(-99.)
    # Infiltration rate of previous day
    RINold = Float(-99)
    # Fraction of non-infiltrating rainfall as function of storm size
    NINFTB = Instance(Afgen)
    # Flag indicating crop present or not
    in_crop_cycle = Bool(False)
    # Flag indicating that a crop was started or finished and therefore the depth
    # of the root zone may have changed, required a redistribution of water
    # between the root zone and the lower zone
    rooted_layer_needs_reset = Bool(False)
    # placeholder for irrigation
    _RIRR = Float(0.)
    # default depth of upper layer (root zone depth)
    DEFAULT_RD = Float(10.)
    # Increments on WLOW due to state updates
    _increments_W = List()

    class Parameters(ParamTemplate):
        # Soil parameters
        SMFCF  = Float(-99.)
        SM0    = Float(-99.)
        SMW    = Float(-99.)
        CRAIRC = Float(-99.)
        SOPE   = Float(-99.)
        KSUB   = Float(-99.)
        RDMSOL = Float(-99.)
        # Site parameters
        IFUNRN = Float(-99.)
        SSMAX  = Float(-99.)
        SSI    = Float(-99.)
        WAV    = Float(-99.)
        NOTINF = Float(-99.)

    class StateVariables(StatesTemplate):
        SM = Float(-99.)
        SS = Float(-99.)
        SSI = Float(-99.)
        W  = Float(-99.)
        WI = Float(-99.)
        WLOW  = Float(-99.)
        WLOWI = Float(-99.)
        WWLOW = Float(-99.)
        # Summation variables 
        WTRAT  = Float(-99.)
        EVST   = Float(-99.)
        EVWT   = Float(-99.)
        TSR    = Float(-99.)
        RAINT  = Float(-99.)
        WDRT   = Float(-99.)
        TOTINF = Float(-99.)
        TOTIRR = Float(-99.)
        PERCT  = Float(-99.)
        LOSST  = Float(-99.)
        # Checksums for rootzone (RT) and total system (TT)
        WBALRT = Float(-99.)
        WBALTT = Float(-99.)
        DSOS = Int(-99)

    class RateVariables(RatesTemplate):
        EVS   = Float(-99.)
        EVW   = Float(-99.)
        WTRA  = Float(-99.)
        RIN   = Float(-99.)
        RIRR  = Float(-99.)
        PERC  = Float(-99.)
        LOSS  = Float(-99.)
        DW    = Float(-99.)
        DWLOW = Float(-99.)
        DTSR = Float(-99.)
        DSS = Float(-99.)
        DRAINT = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: ParameterProvider containing all parameters
        """

        # Check validity of maximum soil moisture amount in topsoil (SMLIM)
        SMLIM = limit(parvalues["SMW"], parvalues["SM0"],  parvalues["SMLIM"])

        if SMLIM != parvalues["SMLIM"]:
            msg = "SMLIM not in valid range, changed from %f to %f."
            self.logger.warn(msg % (parvalues["SMLIM"], SMLIM))

        # Assign parameter values            
        self.params = self.Parameters(parvalues)
        p = self.params
        
        # set default RD to 10 cm, also derive maximum depth and old rooting depth
        RD = self.DEFAULT_RD
        RDM = max(RD, p.RDMSOL)
        self.RDold = RD
        self.RDM = RDM
        
        # Initial surface storage
        SS = p.SSI
        
        # Initial soil moisture content and amount of water in rooted zone,
        # limited by SMLIM. Save initial value (WI)
        SM = limit(p.SMW, SMLIM, (p.SMW + p.WAV/RD))
        W = SM * RD
        WI = W
        
        # initial amount of soil moisture between current root zone and maximum
        # rootable depth (WLOW). Save initial value (WLOWI)
        WLOW  = limit(0., p.SM0*(RDM - RD), (p.WAV + RDM*p.SMW - W))
        WLOWI = WLOW
        
        # Total water depth in soil column (root zone + subsoil)
        WWLOW = W + WLOW

        # soil evaporation, days since last rain (DLSR) set to 1 if the
        # soil is wetter then halfway between SMW and SMFCF, else DSLR=5.
        self.DSLR = 1. if (SM >= (p.SMW + 0.5*(p.SMFCF-p.SMW))) else 5.

        # Initialize some remaining helper variables
        self.RINold = 0.
        self.in_crop_cycle = False
        self.NINFTB = Afgen([0.0,0.0, 0.5,0.0, 1.5,1.0])

        # Initialize model state variables.       
        self.states = self.StateVariables(kiosk, publish=["SM", "DSOS"], SM=SM, SS=SS,
                           SSI=p.SSI, W=W, WI=WI, WLOW=WLOW, WLOWI=WLOWI,
                           WWLOW=WWLOW, WTRAT=0., EVST=0., EVWT=0., TSR=0.,
                           RAINT=0., WDRT=0., TOTINF=0., TOTIRR=0., DSOS=0,
                           PERCT=0., LOSST=0., WBALRT=-999., WBALTT=-999.)
        self.rates = self.RateVariables(kiosk, publish="EVS")
        
        # Connect to CROP_START/CROP_FINISH signals for water balance to
        # search for crop transpiration values
        self._connect_signal(self._on_CROP_START, signals.crop_start)
        self._connect_signal(self._on_CROP_FINISH, signals.crop_finish)
        # signal for irrigation
        self._connect_signal(self._on_IRRIGATE, signals.irrigate)

        self._increments_W = []

    @prepare_rates
    def calc_rates(self, day, drv):
        s = self.states
        p = self.params
        r = self.rates
        k = self.kiosk

        # Rate of irrigation (RIRR)
        r.RIRR = self._RIRR
        self._RIRR = 0.

        # Transpiration and maximum soil and surface water evaporation rates
        # are calculated by the crop evapotranspiration module.
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if "TRA" not in self.kiosk:
            r.WTRA = 0.
            EVWMX = drv.E0
            EVSMX = drv.ES0
        else:
            r.WTRA = k.TRA
            EVWMX = k.EVWMX
            EVSMX = k.EVSMX

        # Actual evaporation rates
        r.EVW = 0.
        r.EVS = 0.
        if s.SS > 1.:
            # If surface storage > 1cm then evaporate from water layer on
            # soil surface
            r.EVW = EVWMX
        else:
            # else assume evaporation from soil surface
            if self.RINold >= 1:
                # If infiltration >= 1cm on previous day assume maximum soil
                # evaporation
                r.EVS = EVSMX
                self.DSLR = 1.
            else:
                # Else soil evaporation is a function days-since-last-rain (DSLR)
                EVSMXT = EVSMX * (sqrt(self.DSLR + 1) - sqrt(self.DSLR))
                r.EVS = min(EVSMX, EVSMXT + self.RINold)
                self.DSLR += 1

        # Potentially infiltrating rainfall
        if p.IFUNRN == 0:
            RINPRE = (1. - p.NOTINF) * drv.RAIN
        else:
            # infiltration is function of storm size (NINFTB)
            RINPRE = (1. - p.NOTINF * self.NINFTB(drv.RAIN)) * drv.RAIN


        # Second stage preliminary infiltration rate (RINPRE)
        # including surface storage and irrigation
        RINPRE = RINPRE + r.RIRR + s.SS
        if s.SS > 0.1:
            # with surface storage, infiltration limited by SOPE
            AVAIL = RINPRE + r.RIRR - r.EVW
            RINPRE = min(p.SOPE, AVAIL)
            
        RD = self._determine_rooting_depth()
        
        # equilibrium amount of soil moisture in rooted zone
        WE = p.SMFCF * RD
        # percolation from rooted zone to subsoil equals amount of
        # excess moisture in rooted zone, not to exceed maximum percolation rate
        # of root zone (SOPE)
        PERC1 = limit(0., p.SOPE, (s.W - WE) - r.WTRA - r.EVS)

        # loss of water at the lower end of the maximum root zone
        # equilibrium amount of soil moisture below rooted zone
        WELOW = p.SMFCF * (self.RDM - RD)
        r.LOSS = limit(0., p.KSUB, (s.WLOW - WELOW + PERC1))

        # percolation not to exceed uptake capacity of subsoil
        PERC2 = ((self.RDM - RD) * p.SM0 - s.WLOW) + r.LOSS
        r.PERC = min(PERC1, PERC2)

        # adjustment of infiltration rate
        r.RIN = min(RINPRE, (p.SM0 - s.SM)*RD + r.WTRA + r.EVS + r.PERC)
        self.RINold = r.RIN

        # rates of change in amounts of moisture W and WLOW
        r.DW = r.RIN - r.WTRA - r.EVS - r.PERC
        r.DWLOW = r.PERC - r.LOSS

        # Check if DW creates a negative value of W
        # If so, reduce EVS to reach W == 0
        Wtmp = s.W + r.DW
        if Wtmp < 0.0:
            r.EVS += Wtmp
            assert r.EVS >= 0., "Negative soil evaporation rate on day %s: %s" % (day, r.EVS)
            r.DW = -s.W

        # Computation of rate of change in surface storage and surface runoff
        # SStmp is the layer of water that cannot infiltrate and that can potentially
        # be stored on the surface. Here we assume that RAIN_NOTINF automatically
        # ends up in the surface storage (and finally runoff).
        SStmp = drv.RAIN + r.RIRR - r.EVW - r.RIN
        # rate of change in surface storage is limited by SSMAX - SS
        r.DSS = min(SStmp, (p.SSMAX - s.SS))
        # Remaining part of SStmp is send to surface runoff
        r.DTSR = SStmp - r.DSS
        # incoming rainfall rate
        r.DRAINT = drv.RAIN

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        p = self.params
        r = self.rates
        
        # INTEGRALS OF THE WATERBALANCE: SUMMATIONS AND STATE VARIABLES

        # total transpiration
        s.WTRAT += r.WTRA * delt

        # total evaporation from surface water layer and/or soil
        s.EVWT += r.EVW * delt
        s.EVST += r.EVS * delt

        # totals for rainfall, irrigation and infiltration
        s.RAINT += r.DRAINT * delt
        s.TOTINF += r.RIN * delt
        s.TOTIRR += r.RIRR * delt

        # Update surface storage and total surface runoff (TSR)
        s.SS += r.DSS * delt
        s.TSR += r.DTSR * delt

        # amount of water in rooted zone
        s.W += r.DW * delt
        assert s.W >= 0., "Negative amount of water in root zone on day %s: %s" % (day, s.W)

        # total percolation and loss of water by deep leaching
        s.PERCT += r.PERC * delt
        s.LOSST += r.LOSS * delt

        # amount of water in unrooted, lower part of rootable zone
        s.WLOW += r.DWLOW * delt
        # total amount of water in the whole rootable zone
        s.WWLOW = s.W + s.WLOW * delt

        # CHANGE OF ROOTZONE SUBSYSTEM BOUNDARY

        # First get the actual rooting depth
        RD = self._determine_rooting_depth()
        RDchange = RD - self.RDold
        self._redistribute_water(RDchange)

        # mean soil moisture content in rooted zone
        s.SM = s.W/RD

        # Accumulate days since oxygen stress, but only if a crop is present
        if s.SM >= (p.SM0 - p.CRAIRC):  # and self.in_crop_cycle:
            s.DSOS += 1
        else:
            s.DSOS = 0

        # save rooting depth
        self.RDold = RD

    @prepare_states
    def finalize(self, day):
        
        s = self.states
        p = self.params

        # Checksums waterbalance for systems without groundwater
        # for rootzone (WBALRT) and whole system (WBALTT)
        # The sum of all increments made are added to ensure a closing waterbalance
        s.WBALRT = s.TOTINF + s.WI + s.WDRT - s.EVST - s.WTRAT - s.PERCT - s.W + sum(self._increments_W)
        s.WBALTT = (s.SSI + s.RAINT + s.TOTIRR + s.WI - s.W + sum(self._increments_W) +
                    s.WLOWI - s.WLOW - s.WTRAT - s.EVWT - s.EVST - s.TSR - s.LOSST - s.SS)

        if abs(s.WBALRT) > 0.0001:
            msg = "Water balance for root zone does not close."
            raise exc.WaterBalanceError(msg)

        if abs(s.WBALTT) > 0.0001:
            msg = "Water balance for complete soil profile does not close.\n"
            msg += ("Total INIT + IN:   %f\n" % (s.WI + s.WLOWI + s.SSI + s.TOTIRR +
                                                 s.RAINT))
            msg += ("Total FINAL + OUT: %f\n" % (s.W + s.WLOW + s.SS + s.EVWT + s.EVST +
                                                 s.WTRAT + s.TSR + s.LOSST))
            raise exc.WaterBalanceError(msg)
        
        # Run finalize on the subSimulationObjects
        SimulationObject.finalize(self, day)
    
    def _determine_rooting_depth(self):
        """Determines appropriate use of the rooting depth (RD)

        This function includes the logic to determine the depth of the upper (rooted)
        layer of the water balance. See the comment in the code for a detailed description.
        """
        if "RD" in self.kiosk:
            return self.kiosk["RD"]
        else:
            # Hold RD at default value
            return self.DEFAULT_RD

    def _redistribute_water(self, RDchange):
        """Redistributes the water between the root zone and the lower zone.

        :param RDchange: Change in root depth [cm] positive for downward growth,
                         negative for upward growth

        Redistribution of water is needed when roots grow during the growing season
        and when the crop is finished and the root zone shifts back from the crop rooted
        depth to the default depth of the upper (rooted) layer of the water balance.
        Or when the initial rooting depth of a crop is different from the default one used
        by the water balance module (10 cm)
        """
        s = self.states
        p = self.params
        
        WDR = 0.
        if RDchange > 0.001:
            # roots grow down by more than 0.001 cm
            # move water from previously unrooted zone and add to new rooted zone
            WDR = s.WLOW * RDchange/(p.RDMSOL - self.RDold)
            # Take minimum of WDR and WLOW to avoid negative WLOW due to rounding
            WDR = min(s.WLOW, WDR)
        else:
            # roots disappear upwards by more than 0.001 cm (especially when crop disappears)
            # move water from previously rooted zone and add to new unrooted zone
            WDR = s.W * RDchange/self.RDold

        if WDR != 0.:
            # reduce amount of water in subsoil
            s.WLOW -= WDR
            # increase amount of water in root zone
            s.W += WDR
            # total water add to rootzone by root zone reset
            s.WDRT += WDR

    def _on_CROP_START(self):
        self.in_crop_cycle = True
        self.rooted_layer_needs_reset = True

    def _on_CROP_FINISH(self):
        self.in_crop_cycle = False
        self.rooted_layer_needs_reset = True

    def _on_IRRIGATE(self, amount, efficiency):
        self._RIRR = amount * efficiency

    def _set_variable_SM(self, nSM):
        """Force the model states based on the given soil moisture value.

        This implies that besides the root-zone soil content, also the
        water available in the root zone (W) must be updated since SM
        is derived from W.

        Further, the increment made to W is added to self._increments_W
        in order to ensure that the water balance still closes.
        """
        s = self.states

        # old values
        oSM = s.SM
        oW = s.W
        # new values
        nW = nSM/oSM * s.W

        # update states
        s.W = nW
        s.SM = nSM
        s.WWLOW = s.WLOW + s.W

        # Store increments on W
        self._increments_W.append(nW - oW)

        # Return increments on all variables
        return {"W": nW - oW, "SM": nSM - oSM}


class WaterbalanceFDSnow(SimulationObject):
    """SimulationObject combining the SnowMAUS and WaterbalanceFD objects.

    Note that the snow module and the water balance are currently treated as
    independent simulations: we are just accumulating a snow pack while
    the soil water balance accumulates the rainfall as if there is no snow.
    This needs to be changed, but the snow module then needs to be integrated
    in a more generic surface storage simulation object that includes all the
    options for storing water on the soil surface: water (ponding), snow and
    also maybe interception by the canopy which is currently lacking.
    """
    waterbalance = Instance(SimulationObject)
    snowcover = Instance(SimulationObject)

    # use observed snow depth
    use_observed_snow_depth = Bool(False)
    # Helper variable for observed snow depth
    _SNOWDEPTH = Float()

    class StateVariables(StatesTemplate):
        SNOWDEPTH = Float()

    def initialize(self, day, kiosk, parvalues):
        self.waterbalance = WaterbalanceFD(day, kiosk, parvalues)

        # Determine of observed or simulated snow depth should be used
        if "ISNOWSRC" not in parvalues:
            msg = "Parameter for selecting observed(0)/simulated(1) snow depth ('ISNOWSRC') missing!"
            raise exc.ParameterError(msg)
        else:
            self.use_observed_snow_depth = True if parvalues["ISNOWSRC"] == 0 else False

        if self.use_observed_snow_depth:
            self.states = self.StateVariables(kiosk, SNOWDEPTH=0.)
        else:
            self.snowcover = SnowMAUS(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.waterbalance.calc_rates(day, drv)
        if self.use_observed_snow_depth:
            self._SNOWDEPTH = drv.SNOWDEPTH
        else:
            self.snowcover.calc_rates(day, drv)

    @prepare_states
    def integrate(self, day, delt=1.0):
        self.waterbalance.integrate(day, delt)
        if self.use_observed_snow_depth:
            self.states.SNOWDEPTH = self._SNOWDEPTH
        else:
            self.snowcover.integrate(day, delt)
