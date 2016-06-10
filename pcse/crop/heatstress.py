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


# -------------------------------------------------------------------------------
class WOFOST_Sink_Dynamics(SimulationObject):
    """Implements the algorithm for heatstress around flowering in WOFOST.
     1/ The number of sinks is established as function of the total leaf and
        stem weight at anthesis.
     2/ The number of sinks is reduced by a factor for high temperatures. This
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
        IHEAT = Float(-99.)  # Switch for heat sensitivity of sink.
        NUMGA = Float(-99.)  # Variable A in relation between TAGP at anthesis anthesis and sink dimension.
        NUMGB = Float(-99.)  # Variable B in relation between TAGP at anthesis anthesis and sink dimension.
        RDGRTB = AfgenTrait()  # High temperature response function for grain formation as a function of day-time temp.

    # -------------------------------------------------------------------------------
    class RateVariables(RatesTemplate):
        DAYTEMP = Float(-99.)  # Increase in day-time temperature sum.

    # -------------------------------------------------------------------------------
    class StateVariables(StatesTemplate):
        CMDTEMP = Float(-99.)  # Sum of the day-time temperature.
        TAGBSF = Float(-99.)  # Total weight of leaves and stems at anthesis.
        NUMGR = Float(-99.)  # number of sinks (per ha).

        # States which register phenological events
        DSB = Instance(datetime.date)  # Day where temp sensitivity begins.
        DSE = Instance(datetime.date)  # Day where temp sensitivity ends.

        HSTAGE = Enum([None, "not sensitive", "temp sensitive", "end temp sensitivity"])

    # ---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):

        HSTAGE = "not sensitive"
        DSB = None
        DSE = None

        self.rates = self.RateVariables(kiosk, publish="DAYTEMP")
        self.params = self.Parameters(parvalues)
        self.states = self.StateVariables(kiosk, publish="NUMGR", NUMGR=0, TAGBSF=0,
                                          CMDTEMP=0., DSB=DSB, DSE=DSE, HSTAGE=HSTAGE)

    # ---------------------------------------------------------------------------

    @prepare_rates
    def calc_rates(self, day, drv):
        # set the daytime temperature as rate variable
        r = self.rates
        r.DAYTEMP = drv.DTEMP

    # ---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day):
        """Updates the state variable and checks for phenologic stages
        """

        p = self.params
        r = self.rates
        s = self.states

        DVS = self.kiosk["DVS"]
        TWST = self.kiosk["TWST"]
        TWLV = self.kiosk["TWLV"]

        number_of_sinks = s.NUMGR

        if DVS >= 1.0 and s.TAGBSF <= 0.:
            # establish total leaf and stem dry weight at anthesis
            s.TAGBSF = TWLV + TWST

        # Calculate the number of sinks as a function of temperature sensitivity
        # and establish if DVS is in the heat sesnsitive period
        if s.HSTAGE == "not sensitive":
            if DVS >= p.DVSHEB:
                self._next_stage(day)

        elif s.HSTAGE == "temp sensitive":
            s.CMDTEMP += r.DAYTEMP
            if DVS >= p.DVSHEF:
                self._next_stage(day)

                # sink reduction factor due to low temperature
                # FTEMP_LOW = p.TMGTB(r.DAYTEMP)

                # number of sinks (per ha) as determined by total leaf and stem
                # dry weight at anthesis
                number_of_sinks = self._number_of_sinks()


        elif s.HSTAGE == "end temp sensitivity":
            dummy = None

        else:  # Problem no heat stage defined
            msg = "No HSTAGE defined in the heat stress around flowering module."
            raise exc.PCSEError(msg)

        s.NUMGR = number_of_sinks

    # ---------------------------------------------------------------------------
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

        else:  # Problem no heat stage defined
            msg = "No HSTAGE defined in temp sensitivity around flowering module."
            raise exc.PCSEError(msg)

        msg = "Changed temp sensivity hstage '%s' to '%s' on %s"
        self.logger.info(msg % (current_HSTAGE, s.HSTAGE, day))

    # ---------------------------------------------------------------------------
    def _number_of_sinks(self):
        # determine number of sinks as defined by heat stress around anthesis and
        # total weight (death+alive) at anthesis
        p = self.params
        s = self.states

        if p.IHEAT == 1:
            delta = s.DSE - s.DSB
            try:
                MDTEMP = s.CMDTEMP / float(delta.days)
            except ZeroDivisionError:
                msg = "Could not calculate heatstress. No days counted"
                raise exc.PCSEError(msg)
            FHSTRESS = p.RDGRTB(MDTEMP)
        else:
            FHSTRESS = 1.0

        return FHSTRESS * FHSTRESS * (p.NUMGA + p.NUMGB * s.TAGBSF)