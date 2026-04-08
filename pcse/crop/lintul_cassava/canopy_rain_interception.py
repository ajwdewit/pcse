# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import SimulationObject, ParamTemplate, RatesTemplate, StatesTemplate, ParameterProvider
from traitlets_pcse import Float

class canopy_rain_interception(SimulationObject):
    """Class to simulate  rain interception by the canopy in the LINTUL Cassava model

    Simulates the daily amount of rain water that is intercepted by the canopy. Depending on the precipitation rate,
    the canopy either intercepts all rain water or it participates a maximum amount per unit of leaf area index.

    **Simulation parameters**

    ==============  ==============================================  ======  ===============================
     Name           Description                                     Type     Unit
    ==============  ==============================================  ======  ===============================
    FRACRNINTC      Maximum daily amount of water that can be
                    intercepted by the canopy per unit of leaf
                    area index                                      SCr     cm water m2 ground m-2 leaf d-1
    ==============  ==============================================  ======  ===============================

    **Rate variables**

    ==============  ==============================================  ======  ===============================
     Name           Description                                     Pbl     Unit
    ==============  ==============================================  ======  ===============================
    RNINTC          Rate of rain interception by the canopy         Y        cm water d-1
    ==============  ==============================================  ======  ===============================

    **Auxillary variables**

    None

    **State variables**

    None
    """

    class Parameters(ParamTemplate):
        FRACRNINTC = Float()

    class RateVariables(RatesTemplate):
        RNINTC = Float()

    class StateVariables(StatesTemplate):
        pass

    def initialize(self, day, kiosk, parvalues):
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk,
                                        publish = ["RNINTC"])
        self.states = self.StateVariables(kiosk,
                                          publish = [])

    def calc_rates(self,  day, drv, delt=1):
        k = self.kiosk
        p = self.params
        r = self.rates

        # Interception of the canopy, depends on the amount of rainfall and the LAI.
        RTRAIN = drv.RAIN / delt  # cm d-1           : rain rate
        RNINTC = min(RTRAIN, (p.FRACRNINTC * k.LAI))  # cm d-1
        r.RNINTC = RNINTC

    def integrate(self, day, drv, delt = 1):
        pass


