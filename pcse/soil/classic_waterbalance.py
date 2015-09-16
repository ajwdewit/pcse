# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Python implementations of the WOFOST waterbalance modules for simulation
of potential production (`WaterbalancePP`) and water-limited production
(`WaterbalanceFD`)under freely draining conditions.
"""
from math import sqrt

from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool
from ..decorators import prepare_rates, prepare_states
from ..util import limit, Afgen, merge_dict
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals
from .. import exceptions as exc
from .snowmaus import SnowMAUS


class WaterbalancePP(SimulationObject):
    """Fake waterbalance for simulation under potential production.
    """
    
    class Parameters(ParamTemplate):
        SMFCF = Float(-99.)

    class StateVariables(StatesTemplate):
        SM = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """    
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param soildata: dictionary with WOFOST soildata key/value pairs
    
        This waterbalance keeps the soil moisture always at field capacity. Therefore   
        `WaterbalancePP` has only one parameter (`SMFCF`: the field capacity of the
        soil) and one state variable (`SM`: the volumetric soil moisture).
        """
        self.params = self.Parameters(parvalues)
        self.states = self.StateVariables(kiosk, SM=self.params.SMFCF,
                                          publish="SM")
    
    @prepare_rates
    def calc_rates(self, day, drv):
        pass
    
    @prepare_states
    def integrate(self, day):
        self.states.SM = self.params.SMFCF
        

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
    soil profile consists of two layers, the actually rooted  soil and the soil
    immediately below the rooted zone until the maximum rooting depth (soil and
    crop dependent). The extension of the root zone from initial rooting depth
    to maximum rooting depth is described in Root_Dynamics class. From the
    moment that the maximum rooting depth is reached the soil profile is
    described as a one layer system. The class WaterbalanceFD is derived
    from WATFD.FOR in WOFOST7.1
    
    **Simulation parameters:** (provide in crop, soil and sitedata dictionary)
    
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    SMFCF     Field capacity of the soil                       SSo     -
    SM0       Porosity of the soil                             SSo     -
    SMW       Wilting point of the soil                        SSo     -
    CRAIRC    Soil critical air content (waterlogging)         SSo     -
    SOPE      maximum percolation rate root zone               SSo    |cmday-1|
    KSUB      maximum percolation rate subsoil                 SSo    |cmday-1|
    K0        hydraulic conductivity of saturated soil         SSo    |cmday-1|
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
    IAIRDU    Switch airducts on (1) or off (0)                SCr     - 
    RDMCR     Maximum rooting depth of the crop                SCr      cm
    RDI       Initial rooting depth of the crop                SCr      cm
    ======== =============================================== =======  ==========

    **State variables:**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    SM        Volumetric moisture content in root zone          Y    -
    SS        Surface storage (layer of water on surface)       N    cm
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
    RAINT     Total amount of rainfall                          N    cm
    WDRT      Amount of water added to root zone by increase    N    cm
              of root growth
    TOTINF    Total amount of infiltration                      N    cm
    TOTIRR    Total amount of irrigation (not implemented       N    cm
              yet)
    PERCT     Total amount of water percolating from rooted     N    cm
              zone to subsoil
    LOSST     Total amount of water lost to deeper soil         N    cm
    WBALRT    Checksum for root zone waterbalance. Will be      N    cm
              calculated within `finalize()`, abs(WBALRT) >
              0.0001 will raise a WaterBalanceError
    WBALTT    Checksum for total waterbalance. Will be          N    cm
              calculated within `finalize()`, abs(WBALTT) >
              0.0001 will raise a WaterBalanceError
    =======  ================================================= ==== ============

    **Rate variables:**

    ======== ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    ======== ================================================= ==== ============
    EVS      Actual evaporation rate from soil                  N    |cmday-1|
    EVW      Actual evaporation rate from water surface         N    |cmday-1|
    WTRA     Actual transpiration rate from plant canopy,       N    |cmday-1|
             is directly derived from the variable "TRA" in
             the evapotranspiration module 
    RAIN     Rainfall rate for current day                      N    |cmday-1|
    RIN      Infiltration rate for current day                  N    |cmday-1|
    RIRR     Irrigation rate for current day. Is always set     N    |cmday-1|
             to zero as irrigation is not handled currently
    PERC     Percolation rate to non-rooted zone                N    |cmday-1|
    LOSS     Rate of water loss to deeper soil                  N    |cmday-1|
    DW       Change in amount of water in rooted zone as a      N    |cmday-1|
             result of infiltration, transpiration and
             evaporation.
    DWLOW    Change in amount of water in subsoil               N    |cmday-1|
    ======== ================================================= ==== ============
    
    
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
    # Flag indicating that a crop was removed and therefore the thickness 
    # of the rootzone shift back to its initial value (params.RDI)
    rooted_layer_needs_reset = Bool(False)
    # placeholder for irrigation
    _RIRR = Float(0.)

    class Parameters(ParamTemplate):
        # Soil parameters
        SMFCF  = Float(-99.)
        SM0    = Float(-99.)
        SMW    = Float(-99.)
        CRAIRC = Float(-99.)
        SOPE   = Float(-99.)
        KSUB   = Float(-99.)
        K0     = Float(-99.)
        RDMSOL = Float(-99.)
        # Site parameters
        IFUNRN = Float(-99.)
        SSMAX  = Float(-99.)
        SSI    = Float(-99.)
        WAV    = Float(-99.)
        NOTINF = Float(-99.)
        # crop parameters
        IAIRDU = Float(-99.)
        RDMCR  = Float(-99.)
        RDI    = Float(-99.)

    class StateVariables(StatesTemplate):
        SM = Float(-99.)
        SS = Float(-99.)
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

    class RateVariables(RatesTemplate):
        EVS   = Float(-99.)
        EVW   = Float(-99.)
        WTRA  = Float(-99.)
        RAIN  = Float(-99.)
        RIN   = Float(-99.)
        RIRR  = Float(-99.)
        PERC  = Float(-99.)
        LOSS  = Float(-99.)
        DW    = Float(-99.)
        DWLOW = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :param soildata: dictionary with WOFOST soildata key/value pairs
        :param sitedata: dictionary with WOFOST sitedata key/value pairs
        """

        # Check validity of maximum soil moisture amount in topsoil (SMLIM)
        SMLIM = parvalues["SMLIM"]
        if parvalues["IAIRDU"] == 1: # applicable only for flooded rice crops
            SMLIM = parvalues["SM0"]
        else:
            SMLIM = limit(parvalues["SMW"], parvalues["SM0"], SMLIM)

        if SMLIM != parvalues["SMLIM"]:
            msg = "SMLIM not in valid range, changed from %f to %f."
            self.logger.warn(msg % (parvalues["SMLIM"], SMLIM))

        # Assign parameter values            
        self.params = self.Parameters(parvalues)
        p = self.params
        
        # Current, maximum and old rooting depth
        RD = p.RDI
        RDM = max(p.RDI, min(p.RDMSOL, p.RDMCR))
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
        self.states = self.StateVariables(kiosk, publish="SM", SM=SM, SS=SS,
                           W=W, WI=WI, WLOW=WLOW, WLOWI=WLOWI, WWLOW=WWLOW,
                           WTRAT=0., EVST=0., EVWT=0., TSR=0.,
                           RAINT=0., WDRT=0., TOTINF=0., TOTIRR=0.,
                           PERCT=0., LOSST=0., WBALRT=-999., WBALTT=-999.)
        self.rates = self.RateVariables(kiosk)
        
        # Connect to CROP_EMERGED/CROP_FINISH signals for water balance to
        # search for crop transpiration values
        self._connect_signal(self._on_CROP_START, signals.crop_start)
        self._connect_signal(self._on_CROP_FINISH, signals.crop_finish)
        # signal for irrigation
        self._connect_signal(self._on_IRRIGATE, signals.irrigate)

    @prepare_rates
    def calc_rates(self, day, drv):
        s = self.states
        p = self.params
        r = self.rates

        # Rate of irrigation (RIRR) set to zero
        r.RIRR = self._RIRR
        self._RIRR = 0.
        
        # Rainfall rate
        r.RAIN = drv.RAIN

        # Transpiration and maximum soil and surface water evaporation rates
        # are calculated by the crop Evapotranspiration module. 
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if "TRA" not in self.kiosk:
            r.WTRA = 0.
            EVWMX = drv.E0
            EVSMX = drv.ES0
        else:
            r.WTRA = self.kiosk["TRA"]
            EVWMX = self.kiosk["EVWMX"]
            EVSMX = self.kiosk["EVSMX"]
        
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
                self.DSLR += 1
                EVSMXT = EVSMX*(sqrt(self.DSLR) - sqrt(self.DSLR-1))
                r.EVS = min(EVSMX, EVSMXT + self.RINold)
        
        # Preliminary infiltration rate (RINPRE)
        if s.SS < 0.1:
            # without surface storage
            if (p.IFUNRN == 0):
                RINPRE = (1. - p.NOTINF)*drv.RAIN + r.RIRR + s.SS

            else:
                RINPRE = (1. - p.NOTINF*self.NINFTB(drv.RAIN)) * drv.RAIN + \
                         r.RIRR + s.SS
        else:
            # with surface storage, infiltration limited by SOPE
            AVAIL = s.SS + (drv.RAIN * (1.-p.NOTINF)) + r.RIRR - r.EVW
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
        LOSS  = limit(0., p.KSUB, (s.WLOW - WELOW + PERC1))
        # for rice water losses are limited to K0/20
        if (p.IAIRDU == 1):
            LOSS = min(LOSS, p.K0/20.)
        r.LOSS = LOSS

        # percolation not to exceed uptake capacity of subsoil
        PERC2 = ((self.RDM -RD) * p.SM0 - s.WLOW) + LOSS
        r.PERC  = min(PERC1,PERC2)

        # adjustment of infiltration rate
        r.RIN = min(RINPRE, (p.SM0 - s.SM)*RD + r.WTRA + r.EVS + r.PERC)
        self.RINold = r.RIN

        # rates of change in amounts of moisture W and WLOW
        r.DW    = r.RIN - r.WTRA - r.EVS - r.PERC
        r.DWLOW = r.PERC - r.LOSS


    @prepare_states
    def integrate(self, day):
        s = self.states
        p = self.params
        r = self.rates
        
        # INTEGRALS OF THE WATERBALANCE: SUMMATIONS AND STATE VARIABLES

        # total transpiration
        s.WTRAT += r.WTRA

        # total evaporation from surface water layer and/or soil
        s.EVWT += r.EVW
        s.EVST += r.EVS

        # totals for rainfall, irrigation and infiltration
        s.RAINT  += r.RAIN
        s.TOTINF += r.RIN
        s.TOTIRR += r.RIRR

        # Update surface storage, any storage > SSMAX goes to total surface
        # runoff (TSR)
        SSPRE = s.SS + (r.RAIN + r.RIRR -r.EVW - r.RIN)
        s.SS  = min(SSPRE, p.SSMAX)
        s.TSR += (SSPRE - s.SS)

        # amount of water in rooted zone
        W_NEW = s.W + r.DW
        if (W_NEW < 0.0):
            # If negative soil water depth, set W to zero and subtract W_NEW
            # from total soil evaporation to keep the balance. 
            # Note: W_NEW is negative here!!
            s.EVST += W_NEW
            s.W = 0.0
        else:
            s.W = W_NEW

        # total percolation and loss of water by deep leaching
        s.PERCT += r.PERC
        s.LOSST += r.LOSS

        # amount of water in unrooted, lower part of rootable zone
        s.WLOW += r.DWLOW
        # total amount of water in the whole rootable zone
        s.WWLOW = s.W + s.WLOW

        # CHANGE OF ROOTZONE SUBSYSTEM BOUNDARY

        # Redefine the rootzone when the crop is finished. As a result there
        # no roots anymore (variable RD) and the rootzone shifts back to its
        # initial depth (params.RDI). As a result the amount of water in the
        # initial rooting depth and the unrooted layer must be redistributed. 
        # Note that his is a rather artificial solution resulting from the fact
        # that the rooting depth is user as to define a layer in the WOFOST
        # water balance.
        if self.rooted_layer_needs_reset is True:
            self._reset_rootzone()

        # calculation of new amount of soil moisture in rootzone by root growth
        RD = self._determine_rooting_depth()            
        if (RD - self.RDold) > 0.001:
            # water added to root zone by root growth, in cm
            WDR = s.WLOW * (RD - self.RDold)/(self.RDM - self.RDold)
            s.WLOW -= WDR

            # total water addition to rootzone by root growth
            s.WDRT += WDR
            # amount of soil moisture in extended rootzone
            s.W += WDR

        # mean soil moisture content in rooted zone
        s.SM = s.W/RD
        # save rooting depth
        self.RDold = RD

    @prepare_states
    def finalize(self, day):
        
        s = self.states
        p = self.params

        # Checksums waterbalance for systems without groundwater
        # for rootzone (WBALRT) and whole system (WBALTT)
        s.WBALRT = s.TOTINF + s.WI + s.WDRT - s.EVST - s.WTRAT - s.PERCT - s.W
        s.WBALTT = (p.SSI + s.RAINT + s.TOTIRR + s.WI - s.W + s.WLOWI - 
                    s.WLOW - s.WTRAT- s.EVWT - s.EVST - s.TSR - s.LOSST - s.SS)
        if abs(s.WBALRT) > 0.0001:
            msg = "Water balance for root zone does not close."
            raise exc.WaterBalanceError(msg)

        if abs(s.WBALTT) > 0.0001:
            msg = "Water balance for complete soil profile does not close.\n"
            msg += ("Total INIT + IN:   %f\n" % (s.WI+s.WLOWI+p.SSI+s.TOTIRR+
                                                 s.RAINT))
            msg += ("Total FINAL + OUT: %f\n" % (s.W+s.WLOW+s.SS+s.EVWT+s.EVST+
                                                 s.WTRAT+s.TSR+s.LOSST))
            raise exc.WaterBalanceError(msg)
        
        # Run finalize on the subSimulationObjects
        SimulationObject.finalize(self, day)
    
    def _determine_rooting_depth(self):
        """Determines appropriate use of the rooting depth (RD)
        """
        p = self.params

        if self.in_crop_cycle is False:

            # Crop finished
            if  "RD" in self.kiosk:
                # Only happens at the final simulation cycle when value for
                # SM still has to be computed. This also applies that a reset
                # of the root zone layer is needed in the next update cycle.
                RD = self.kiosk["RD"]
                self.rooted_layer_needs_reset = True
            else:
                # not in crop cycle hold RD at initial value RDI
                RD = p.RDI
            
        else: # In cropping season
            RD = self.kiosk["RD"]
            
        return RD
    
    def _reset_rootzone(self):
        s = self.states
        p = self.params
        
        self.rooted_layer_needs_reset = False
        
        # water added to the subsoil by root zone reset
        WDR = s.W * (self.RDold - p.RDI)/(self.RDold)
        s.WLOW += WDR

        # total water subtracted from, rootzone by root zone reset
        s.WDRT -= WDR
        # amount of soil moisture in new resetted rootzone
        s.W -= WDR
                
        
    def _on_CROP_START(self):
        self.in_crop_cycle = True
        
    def _on_CROP_FINISH(self):
        self.in_crop_cycle = False

    def _on_IRRIGATE(self, amount, efficiency):
        self._RIRR = amount * efficiency


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

    def initialize(self, day, kiosk, parvalues):
        self.waterbalance = WaterbalanceFD(day, kiosk, parvalues)
        self.snowcover = SnowMAUS(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.waterbalance.calc_rates(day, drv)
        self.snowcover.calc_rates(day, drv)

    def integrate(self, day):
        self.waterbalance.integrate(day)
        self.snowcover.integrate(day)
