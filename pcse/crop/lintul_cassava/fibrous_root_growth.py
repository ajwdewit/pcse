# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float

class fibrous_root_growth(SimulationObject):
    """
    Class to simulate the rooting depth.

    The class calculates the increase in rooting depth from emergence until a maximum rooting depth is reached.

    **Simulation parameters**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Type     Unit
    =================  ==============================================  ======  ===========================
    RDI                Initial rooting depth                           SCr     cm root
    RDMSOL             Maximum rooting depth                           SCr     cm root
    RRDMAX             Maximum growth rate of rooting depth            SCr     cm root
    SMW                Soil moisture content at permanent wilting
                       point                                           SCr     cm root
    =================  ==============================================  ======  ===========================

    **State variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl      Unit
    =================  ==============================================  ======  ===========================
    RD                 Rooting depth                                   Y       cm root
    =================  ==============================================  ======  ===========================

    **Rate variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl      Unit
    =================  ==============================================  ======  ===========================
    RRD                Rate of change of rooting depth                 N        cm root
    =================  ==============================================  ======  ===========================

    """

    class Parameters(ParamTemplate):
        RDI = Float()
        RDMSOL = Float()
        RRDMAX = Float()
        SMW = Float()

    class RateVariables(RatesTemplate):
        RRD = Float()

    class StateVariables(StatesTemplate):
        RD = Float()

    def initialize(self, day, kiosk, parvalues):
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        k = self.kiosk
        p = self.params
        RD = p.RDI
        self.rates = self.RateVariables(kiosk,
                                        publish = ["RRD"])
        self.states = self.StateVariables(kiosk,
                                          publish = ["RD"],
                                          RD = RD)

    def calc_rates(self,  day, drv, delt=1):
        # If the soil water content drops to, or below, wilting point fibrous root growth stops.
        # Root growth continues till the maximum rooting depth is reached.
        # The rooting depth (cm) is calculated from a maximum rate of change in rooting depth,
        # the emergence of the crop and the constraints mentioned above.

        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        if (s.RD-p.RDMSOL < 0) & (k.SM-p.SMW >= 0):
            RRD = p.RRDMAX * k.EMERG  # cm d-1
        else:
            RRD = 0

        r.RRD = RRD

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states
        s.RD += r.RRD * delt
