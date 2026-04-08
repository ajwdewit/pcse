# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float

class dormancy_and_recovery(SimulationObject):
    """
    Class to simulate dormancy and recovery in LINTUL Cassava.

    Simulates the initiation of dormancy, which occurs once the soil water content is lower than the soil water content
    at severe drought and the LAI is low, and the first branching has not yet taken place. It also simulates the end of
    dormancy and the beginning of the recovery phase. Recovery is initiated once the soil water content is above the
    soil water content at severe drought. During this period, storage organ dry matter can be redistributed to the
    leaves. At the end of the recovery period, redistribution stops.

    **Simulation parameters**

    =================  ================================================  ======  ===========================
    Name               Description                                       Type     Unit
    =================  ================================================  ======  ===========================
    DELREDIST          Delay for redistribution of dry matter             SCr      |C| d
    LAI_MIN            Leaf are index below which the dormancy phase
                       can be entered.                                    SCr      |C| d
    RECOV              Fraction of soil moisture content above which
                       the crop recovers from dormancy to the
                       critical soil moisture content                     SCr      cm water cm-1 water
    RREDISTSO          Relative rate of dry matter redistribution
                       from the storage organs to leaves
    SMW                Soil moisture content at wilting point             SCr      cm3 water cm-3 soil
    SO2LV              Conversion rate of storage organ dry matter
                       to leaf dry matter                                 SCr      g DM g-1 DM
    TSUMSBR            Temperature sum at which first branching takes
                       place                                              SCr      |C| d
    TSUMREDISTMAX      Temperature sum of duration of dry matter
                       redistribution                                     SCr      |C| d
    WSOREDISTFRACMAX   Maximum fraction of storage organ dry matter
                       that can be redistributed to the leaves.           SCr      |C| d
    WLVGNEW            Minimum amount of dry matter of new leaves that    
                       can be produced in the redistribution phase.       SCr      g DM m-2 ground
    =================  ================================================  ======  ===========================

    **State variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    DORMTIME           Number of days that the crop was in dormancy    N       d
    DORMTSUM           Temperature sum of dormancy period              N       |C| d
    PUSHREDISTENDTSUM
    PUSHDORMRECTSUM
    PUSHREDISTSUM
    REDISTLVG          Amount of green dry weight that was produced
                       by dry weight redistribution from storage
                       organs to green leaves                          N       g DM m-2 ground
    REDISTSO           Amount of storage organ dry weight that was     N
                       lost due to dry weight redistribution
                       from storage organs to green leaves.                    g DM m-2 ground
    =================  ==============================================  ======  ===========================

    **Rate variables**

    ===================  ==============================================  ======  ===========================
    Name                 Description                                     Pbl     Unit
    ===================  ==============================================  ======  ===========================
    RDORMTIME            Rate at which the number of days in dormancy    
                         changes                                          N        d d-1
    RDORMTSUM            Rate of change of temperature sum for the
                         cperiod in which the rop is in dormancy          N        (|C| d) d-1
    RPUSHDORMRECTSUM                                                      N        -
    RPUSHREDISTSUM                                                        N        -
    RPUSHREDISTENDTSUM                                                    N        -
    RREDISTLVG           Rate of change of green leaf dry weight due      N        g DM m-2 ground d-1
                         to redistribution from storage organs to         Y
                         green leaves                                     Y        g DM m-2 ground d-1
    RREDISTSO            Rate of change of storage organdry weight
                         due to redistribution from storage organs to
                         green leaves
    ===================  ==============================================  ======  ===========================

    **Auxillary variables**
    
    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================    
    DORMANCY           Indicates whether (1) or not (0) the crop is          
                       in the dormancy phase                           Y        -
    PUSHREDIST                                                         Y     
    =================  ==============================================  ======  ===========================    

    """

    class Parameters(ParamTemplate):
        DELREDIST = Float()
        LAI_MIN = Float()
        RECOV = Float()
        RRREDISTSO = Float()
        SO2LV = Float()
        TSUMSBR = Float()
        TSUMREDISTMAX = Float()
        WSOREDISTFRACMAX = Float()
        WLVGNEWN = Float()
        SMW = Float()

    class RateVariables(RatesTemplate):
        DORMANCY = Float()
        RDORMTSUM = Float()
        RPUSHDORMRECTSUM = Float()
        RPUSHREDISTENDTSUM = Float()
        RDORMTIME = Float()
        RREDISTLVG = Float()
        RREDISTSO = Float()
        RPUSHREDISTSUM = Float()
        PUSHREDIST = Float()

    class StateVariables(StatesTemplate):
        DORMTSUM = Float()
        PUSHREDISTENDTSUM = Float()
        PUSHDORMRECTSUM = Float()
        DORMTIME = Float()
        REDISTLVG = Float()
        REDISTSO = Float()
        PUSHREDISTSUM = Float()

    def initialize(self, day, kiosk, parvalues):
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk,
                                        publish = ["DORMANCY",
                                                   "PUSHREDIST",
                                                   "RREDISTLVG",
                                                   "RREDISTSO"])

        DORMTSUM = 0.
        PUSHDORMRECTSUM = 0.
        PUSHREDISTENDTSUM = 0.
        DORMTIME = 0.
        REDISTLVG = 0.
        REDISTSO = 0.
        PUSHREDISTSUM = 0.
        self.states = self.StateVariables(kiosk,
                                          publish = [],
                                          DORMTSUM=DORMTSUM,
                                          PUSHDORMRECTSUM=PUSHDORMRECTSUM,
                                          PUSHREDISTENDTSUM=PUSHREDISTENDTSUM,
                                          DORMTIME=DORMTIME,
                                          REDISTLVG=REDISTLVG,
                                          REDISTSO=REDISTSO,
                                          PUSHREDISTSUM=PUSHREDISTSUM
                                          )

    def calc_rates(self, day, drv, delt=1):
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        # The crop enters the dormancy phase as the soil water content is lower than the soil water content at
        # severe drought and as the LAI is lower than the minimal LAI.
        if (k.SM-k.WCSD <= 0)  & (k.LAI - p.LAI_MIN <= 0):
            dormancy = 1
        else:
            dormancy = 0

        # # The crop goes out of dormancy if the water content is higher than a certain recovery water content and as the
        # # water content is larger than the wilting point soil moisture content.
        if (k.SM - p.RECOV * k.WCCR >= 0) & (k.SM - p.SMW >= 0):
            pushdor = 1
        else:
            pushdor = 0

        # # The redistributed fraction of storage root DM to the leaves.
        if k.WSO == 0:
            WSOREDISTFRAC = 1
        else:
            WSOREDISTFRAC = s.REDISTSO / k.WSO

        # Three push functions are used to determine the redistribution and recovery from dormancy, a final function DORMANCY
        # is used to indicate if the crop is still in dormancy:
        # (1) PUSHREDISTEND: The activation of the PUSHREDISTEND function ends the redistribution phase. Redistribution stops
        # when the redistributed fraction reached the maximum redistributed fraction or when the minimum amount of new leaves
        # is produced after dormancy or when the Tsum during the recovery exceeds the maximum redistribution temperature sum.
        # (2) PUSHREDIST: The activation of the PUSHREDIST function ends the dormancy phase including the delay temperature
        # sum needed for the redistribution of DM.
        # (3) PUSHDORMREC: Indicates if the the crop is still in dormancy. Dormancy can only when the temperature sum of the
        # crop exceeds the temperature sum of the branching.
        if WSOREDISTFRAC - p.WSOREDISTFRACMAX >= 0:
            PUSHREDISTEND1 = 1
        else:
            PUSHREDISTEND1 = 0

        if s.REDISTLVG - p.WLVGNEWN  >= 0:
            PUSHREDISTEND2 = 1
        else:
            PUSHREDISTEND2 = 0

        if s.PUSHREDISTSUM - p.TSUMREDISTMAX >= 0:
            PUSHREDISTEND3 = 1
        else:
            PUSHREDISTEND3 = 0

        if -s.PUSHREDISTSUM >= 0:
            PUSHREDISTEND4 = 0
        else:
            PUSHREDISTEND4 = 1

        PUSHREDISTEND = max(max([PUSHREDISTEND1, PUSHREDISTEND2]), PUSHREDISTEND3 * PUSHREDISTEND4)

        if s.PUSHDORMRECTSUM - p.DELREDIST >= 0:
            PUSHREDIST1 = 1
        else:
            PUSHREDIST1 = 0

        PUSHREDIST2 = (1- PUSHREDISTEND)
        PUSHREDIST = PUSHREDIST1 * PUSHREDIST2

        if -s.DORMTSUM >= 0:
            PUSHDORMREC1 = 0
        else:
            PUSHDORMREC1 = 1

        if k.TSUMCROP - p.TSUMSBR >= 0:
            PUSHDORMREC2 = 1
        else:
            PUSHDORMREC2 = 0

        PUSHDORMREC = pushdor * PUSHDORMREC1 * (1 - PUSHREDIST) * PUSHDORMREC2

        if k.TSUMCROP - p.TSUMSBR >= 0:
            DORMANCY1 = 1
        else:
            DORMANCY1 = 0

        DORMANCY = max(dormancy, PUSHDORMREC) * (1 - PUSHREDIST) * DORMANCY1 # (-)

        # The temperature sums related to the dormancy and recovery periods.
        RDORMTSUM = k.DTEFF * DORMANCY - (s.DORMTSUM / delt) * PUSHREDIST  # Deg. C
        RPUSHDORMRECTSUM = k.DTEFF * PUSHDORMREC - (s.PUSHDORMRECTSUM / delt) * (1 - PUSHDORMREC) * (1 - PUSHREDIST)  # Deg. C
        RPUSHREDISTSUM = k.DTEFF * PUSHREDIST - (s.PUSHREDISTSUM / delt) * PUSHREDISTEND  # Deg. C
        RPUSHREDISTENDTSUM = k.DTEFF * PUSHREDIST - (s.PUSHREDISTENDTSUM / delt) * (1 - PUSHREDISTEND)  # Deg. C

        # No. of days in dormancy
        RDORMTIME = DORMANCY  # d

        if -s.DORMTSUM >= 0:
            RREDISTSO1 = 0
        else:
            RREDISTSO1 = 1

        # Dry matter redistribution after dormancy. The rate of redistribution of the storage roots dry matter to
        # leaf dry matter. A certain fraction is lost for the conversion of storage organs dry matter to leaf dry
        # matter.
        RREDISTSO = p.RRREDISTSO * k.WSO * PUSHREDIST - (s.REDISTSO / delt) * RREDISTSO1  # g DM m-2 d-1
        RREDISTLVG = p.SO2LV * RREDISTSO * (1- DORMANCY)  # g DM m-2 d-1
        RREDISTMAINTLOSS = (1 - p.SO2LV) * RREDISTSO  # g DM m-2 d-1

        r.RDORMTSUM = RDORMTSUM
        r.RPUSHDORMRECTSUM = RPUSHDORMRECTSUM
        r.RPUSHREDISTENDTSUM = RPUSHREDISTENDTSUM
        r.RDORMTIME = RDORMTIME
        r.RREDISTLVG = RREDISTLVG
        r.RREDISTSO = RREDISTSO
        r.RPUSHREDISTSUM = RPUSHREDISTSUM

        r.DORMANCY = DORMANCY
        r.PUSHREDIST = PUSHREDIST

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states
        s.DORMTSUM += delt * r.RDORMTSUM
        s.PUSHDORMRECTSUM += delt * r.RPUSHDORMRECTSUM
        s.PUSHREDISTENDTSUM += delt * r.RPUSHREDISTENDTSUM
        s.DORMTIME += delt * r.RDORMTIME
        s.REDISTLVG += delt * r.RREDISTLVG
        s.REDISTSO += delt * r.RREDISTSO
        s.PUSHREDISTSUM += delt * r.RPUSHREDISTSUM
        pass