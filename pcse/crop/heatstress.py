# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl )and Iwan Supit
# (iwan.supit@wur.nl), June 2016
# Approach based on LINTUL4-vsht made by Joost Wolf

"""Implementation of a models for heatstress around flowering in WOFOST

Classes defined here:
- WOFOST_Sink_Dynamics:
"""
import datetime

from ..traitlets import Float, Int, Instance, Enum, Bool, AfgenTrait
from ..decorators import prepare_rates, prepare_states

from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject, VariableKiosk
from .. import exceptions as exc

#-------------------------------------------------------------------------------
class WOFOST_Sink_Dynamics(SimulationObject):
    """Implements the algorithm for heatstress around flowering in WOFOST.
     1/ The number of sinks is established as function of the total leaf and
        stem weight at anthesis.
     2/ This number of sinks is reduced by a factor for high temperatures. This
        factor depends on the average day time temperature during the sensitive
        period.
     3/ The number of sinks is also reduced by low temperature at the end of the
        sensitive period. This reduction factor depends on the day time temperature
        at DVSHEF.


    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    RDGRTB   Reduction factor of grain formation due        TCr         -
             to heat stress in the period around anthesis   TCr         -
    TMGTB    Reduction factor of grain formation  as
             function. Values from SWHEAT model


    DVSHEB   Beginning of the period with possible temp     SCr        |day|
             sensitivity effects
    DVSHEF   End of the period with possible heat temp    SCr        |day|
             sensitivity effects
    NUMGA    Variable A in relation between TAGP at         SCr        |ha-1|
             anthesis and sink dimension
    NUMGB    Variable B in relation between TAGP at         SCr        |ha-1|
             anthesis and sink dimension
    =======  ============================================= =======  ============

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    CMDTEMP  Temperature sum day-time temperature               N    |C| day
             between DVSHEB and DVSHEF
    DSB      Starting day of heat sensitivity                   N    -
    DSE      End day of heat sensitivity                        N    -
    HSTAGE   Current phenological stage, can take the           N    -
             folowing values:
             `not sensitive|temp sensitive|end temp sensitivty`
    NUMGR    Number of grains                                   Y    -
    TAGBSF   Total leaf and stem weight (dead/alive)           N    |kg ha-1|
             at anthesis
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    DAYTEMP  Daily increase of day-time temperature             N       |C|
    =======  ================================================= ==== ============
    
    **External dependencies:**


    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    TWLV     Total weight of leaves              leaf_dynamics       |kg ha-1|
    TWST     Total weight of stems               stem_dynamics       |kg ha-1|
    =======  =================================== =================  ============

    """

    class Parameters(ParamTemplate):
        DVSHEB = Float(-99.)  # beginning of the period with possible heat stress effects.
        DVSHEF = Float(-99.)  # finalization of the period with possible heat stress effects.
        IHEAT  = Float(-99.)  # Switch for heat sensitivity of sink.
        NUMGA  = Float(-99.)  # Variable A in relation between TAGP at anthesis anthesis and sink dimension.
        NUMGB  = Float(-99.)  # Variable B in relation between TAGP at anthesis anthesis and sink dimension.
        RDGRTB = AfgenTrait()  # High temperature response function for grain formation as a function of day-time temp.
        TMGTB  = AfgenTrait()  # Low temperature response function for grain formation as a function of day-time temp.
    #-------------------------------------------------------------------------------
    class RateVariables(RatesTemplate):
        DAYTEMP = Float(-99.)  # Increase in day-time temperature sum.


    #-------------------------------------------------------------------------------
    class StateVariables(StatesTemplate):
        CMDTEMP   = Float(-99.)  # Sum of the day-time temperature.
        TAGBSF    = Float(-99.)  # Total weight of leaves and stems at anthesis.
        NUMGR     = Float(-99.)  # number of sinks (per ha).

        # States which register phenological events
        DSB = Instance(datetime.date) # Day where temp sensitivity begins.
        DSE = Instance(datetime.date) # Day where temp sensitivity ends.

        HSTAGE = Enum([None, "not sensitive", "temp sensitive", "end temp sensitivity"])

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):

        HSTAGE = "not sensitive"
        DSB = None
        DSE = None
        TAGBSF = NUMGR = 0.

        self.rates  = self.RateVariables(kiosk,publish="DAYTEMP")
        self.params = self.Parameters(parvalues)
        self.states = self.StateVariables(kiosk, publish="NUMGR",NUMGR=NUMGR,TAGBSF=TAGBSF,
                                          CMDTEMP=0.,DSB=DSB, DSE=DSE,HSTAGE=HSTAGE)

    #---------------------------------------------------------------------------

    @prepare_rates
    def calc_rates(self, day, drv):
        # set the daytime temperature as rate variable
        r = self.rates
        r.DAYTEMP = drv.DTEMP


    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day):
        """Updates the state variable and checks for phenologic stages
        """

        p = self.params
        r = self.rates
        s = self.states

        DVS  = self.kiosk["DVS"]
        TWST = self.kiosk["TWST"]
        TWLV = self.kiosk["TWLV"]


        # Calculate the number of sinks as a function of temperature sensitivity
        if s.HSTAGE == "not sensitive":
            s.CMDTEMP = 0.
            s.NUMGR   = 0.
            if DVS >= p.DVSHEB:
                self._next_stage(day)

            if DVS >= 1.0 and s.TAGBSF is None:
                # establish total leaf and stem dry weight at anthesis
                s.TAGBSF = TWLV + TWST

        elif s.HSTAGE == "temp sensitive":
            s.CMDTEMP += r.DAYTEMP
            s.NUMGR = 0.
            if DVS >= p.DVSHEF:
                self._next_stage(day)

                # sink reduction factor due to low temperature
                TMG = p.TMGTB(r.DAYTEMP)

                # number of sinks (per ha) as determined by total leaf and stem
                # dry weight at anthesis
                if p.IHEAT == 1:
                    s.NUMGR = self._number_of_sinks_heatstress(TMG)
                else:
                    s.NUMGR = self._number_of_sinks_no_heatstress(TMG)
        elif s.HSTAGE == "end temp sensitivity":
            s.NUMGR=s.NUMGR



        else: # Problem no heat stage defined
            msg = "No HSTAGE defined in the heat stress around flowering module."
            raise exc.PCSEError(msg)


    #---------------------------------------------------------------------------
    def _next_stage(self, day):
        """Moves states.STAGE to the next temp sensitivity stage"""
        s = self.states

        current_HSTAGE = s.HSTAGE
        if s.HSTAGE == "not sensitive":
            s.HSTAGE = "temp sensitive"
            s.DSB = day
            
        elif s.HSTAGE == "temp sensitive":
            s.HSTAGE = "end temp sensitivity"
            s.DSE = day

        else: # Problem no heat stage defined
            msg = "No HSTAGE defined in temp sensitivity around flowering module."
            raise exc.PCSEError(msg)
        
        msg = "Changed temp sensivity hstage '%s' to '%s' on %s"
        self.logger.info(msg % (current_HSTAGE, s.HSTAGE, day))

    #---------------------------------------------------------------------------
    def _number_of_sinks_heatstress(self, TMG):
        # determine number of sinks as defined by heat stress around anthesis and
        # total weight (death+alive) at anthesis
        p = self.params
        s = self.states

        delta  = s.DSE - s.DSB
        try:
            MDTEMP = s.CMDTEMP / float(delta.days)
        except ZeroDivisionError:
            msg = "Could not calculate heatstress. No days counted"
            raise exc.PCSEError(msg)

        return p.RDGRTB(MDTEMP) * TMG * (p.NUMGA + p.NUMGB * s.TAGBSF)

    def _number_of_sinks_no_heatstress(self, TMG):
        # determine number of sinks without heat stress around anthesis
        p = self.params
        s = self.states

        # Note that suboptimal temperatures at the end of the heat sensitive period
        # may also reduce number of sinks
        return TMG * (p.NUMGA + p.NUMGB * s.TAGBSF)



