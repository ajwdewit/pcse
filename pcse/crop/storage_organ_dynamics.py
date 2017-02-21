# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
# Adapted for translocation, sink limitation and
# heat stress by Iwan Supit (iwan.supit@wur.nl), June 2016

from ..traitlets import Float, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk


class WOFOST_Storage_Organ_Dynamics(SimulationObject):
    """Implementation of storage organ dynamics.

    Storage organs are the most simple component of the plant in WOFOST and
    consist of a static pool of biomass. Growth of the storage organs is the
    result of assimilate partitioning. Death of storage organs is not
    implemented and the corresponding rate variable (DRSO) is always set to
    zero.

    Pods are green elements of the plant canopy and can as such contribute
    to the total photosynthetic active area. This is expressed as the Pod
    Area Index which is obtained by multiplying pod biomass with a fixed
    Specific Pod Area (SPA).

    If the storage organs growth also depends on sink limitation (ISINK==1),
    low temperature reduces storage organs growth. In case the
    sink_limitation is less than the soure_limitation the translocation
    amount and consequently also the the net increase of the net stem
    biomass should be adapted.

    **Simulation parameters**

    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    TDWI     Initial total crop dry weight                  SCr      |kg ha-1|
    SPA      Specific Pod Area                              SCr      |ha kg-1|
    GRMX     Maximal grain mass one grain                   SCr      |kg|
    PGRIG    Potential grain formation per grain            SCr      |kg|
    TMGTB    Reduction factor of grain formation  as        TCr       -
             function. Values from SWHEAT model
    =======  ============================================= =======  ============

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    PAI      Pod Area Index                                     Y     -
    WSO      Weight of living storage organs                    Y     |kg ha-1|
    DWSO     Weight of dead storage organs                      N     |kg ha-1|
    TWSO     Total weight of storage organs                     Y     |kg ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    GRSO     Growth rate storage organs                         N   |kg ha-1 d-1|
    DRSO     Death rate storage organs                          N   |kg ha-1 d-1|
    GWSO     Net change in storage organ biomass                N   |kg ha-1 d-1|
    TRANSL   Weight that is translocated to storage organs      N   |kg ha-1 d-1|
    REDUCTL  Reduction of TRANSL resulting from sink            Y   |kg ha-1 d-1|
             limitation
    =======  ================================================= ==== ============

    **Signals send or handled**

    None

    **External dependencies**

    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    ADMI     Above-ground dry matter             CropSimulation     |kg ha-1 d-1|
             increase
    FO       Fraction biomass to storage organs  DVS_Partitioning    -
    FR       Fraction biomass to roots           DVS_Partitioning    -
    TRANSL   Weight that is translocated to      stem_dynamics      |kg ha-1 d-1|
             storage organs
    NUMGR    Number of sinks (grains, etc.)      heatstress          -
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):
        SPA = Float(-99.)
        TDWI = Float(-99.)
        ISINK = Float(-99.)
        PGRIG = Float(-99.)
        GRMX = Float(-99.)
        TMGTB = AfgenTrait()  # Low temperature response function for grain formation as a function of day-time temp.

    class StateVariables(StatesTemplate):
        WSO = Float(-99.)  # Weight living storage organs
        DWSO = Float(-99.)  # Weight dead storage organs
        TWSO = Float(-99.)  # Total weight storage organs
        PAI = Float(-99.)  # Pod Area Index

    class RateVariables(RatesTemplate):
        GRSO = Float(-99.)
        DRSO = Float(-99.)
        GWSO = Float(-99.)
        REDUCTL = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk,publish=["REDUCTL"])
        self.kiosk = kiosk

        # INITIAL STATES
        params = self.params
        # Initial storage organ biomass
        FO = self.kiosk["FO"]
        FR = self.kiosk["FR"]
        WSO = (params.TDWI * (1 - FR)) * FO
        DWSO = 0.
        TWSO = WSO + DWSO
        # Initial Pod Area Index
        PAI = WSO * params.SPA

        self.states = self.StateVariables(kiosk, publish=["TWSO", "WSO", "PAI"],
                                          WSO=WSO, DWSO=DWSO, TWSO=TWSO,
                                          PAI=PAI)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        states = self.states
        params = self.params

        FO = self.kiosk["FO"]
        ADMI = self.kiosk["ADMI"]
        TRANSL = self.kiosk["TRANSL"]

        rates.GRSO = ADMI * FO
        rates.DRSO = 0.0
        rates.GWSO = rates.GRSO - rates.DRSO + TRANSL

        if params.ISINK == 1:
            # set the source growth rate of the storage organs
            GWSO_SR = rates.GWSO
            NUMGR = self.kiosk["NUMGR"]
            rates.REDUCTL = 0.
            if NUMGR > 0.:
                # Low temperature reduction factor for the
                # potential grain growth PGRIG
                FTLOW = params.TMGTB(drv.DTEMP)
                # potential growth of the grains GWSO_PT
                GWSO_PT = NUMGR * FTLOW * params.PGRIG
                # maximal grain mass growth GWSO_MX
                GWSO_MX = max(0.0, NUMGR * params.GRMX - states.WSO)
                # grain  growth limited by both maximal grain mass GWSO_MX and by
                # potential growth of the grains GWSO_PT
                GWSO_SK = min(GWSO_MX, GWSO_PT)
                # source or sink limitation
                rates.GWSO = min(GWSO_SR, GWSO_SK)
                # in case of sink limitation
                if GWSO_SK < GWSO_SR:
                    # adapt translocation rate, not all
                    # assimilates can be stored in storage organs
                    # consequently stem growth should also
                    # be increased with the translocation reduction REDUCTL.
                    rates.REDUCTL = (GWSO_SR - GWSO_SK)


    @prepare_states
    def integrate(self, day, delt=1.0):
        params = self.params
        rates = self.rates
        states = self.states

        # Stem biomass (living, dead, total)
        states.WSO += rates.GWSO
        states.DWSO += rates.DRSO
        states.TWSO = states.WSO + states.DWSO

        # Calculate Pod Area Index (SAI)
        states.PAI = states.WSO * params.SPA






