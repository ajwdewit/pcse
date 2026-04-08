# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float
from pcse.util import AfgenTrait

class leaf_senescence(SimulationObject):
    """
    Class to simulate the growth of green leaf area index.

    This class calculates the various mechanisms that lead to leaf senescence and, because of that, a reduction
    of the leaf area index. Green leaf tissue can either die due to 1) Ageing, 2) Leaf shedding in case of
    too much self shading, 3) severe drought, or 4) nutrient limitation.

    **Simulation parameters**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Type     Unit
    =================  ==============================================  ======  ===========================
    FASTRANSLO         Fraction of senescened leaf dry matterr that
                       is transferred to the storage organs before
                       shedding                                        SCr     g DM g-1 DM
    FRACTLLFENHSH      Fraction of maximum leaf age above which leaf
                       shedding is enhanced.                           SCr      -
    FRACSLATB          Fraction of actual specific leaf area to the
                       maximum specific leaf area as a function of
                       temperature sum                                 TCr      -
    LAICR              Critical leaf area index above which leaf
                       shedding is induced                             SCr      m2 leaf m-2 ground
    RDRB               Relative death rates of leaves in absence
                       of water and nutrient stress                    SCr      d-1
    RDNS               Addition factor for leaf death due to nutrient
                       stress                                          SCr      d-1
    RDRT               Relative death rate of leaves due to
                       temperature                                     SCr      d-1
    RDRSHM             Relative death rate of leaves due to self-
                       shading                                         SCr      d-1
    SLA_MAX            Maximum specific leaf area                      SCr      m2 leaf kg-1 leaf
    TSUMLIFE           Temperature sum of leaf life above which
                       leaf senescence due to ageing starts            SCr      |C| d
    WCWET              Soil moisture content above which oxyen stress
                       occurs.                                         SCr      cm3 water cm-3 ground
    =================  ==============================================  ======  ===========================

    **State variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    TSUMCROPLEAFAGE    Physiological leaf age                          N       |C| d
    WLVD               Dry weight of senescened leaf dry matter        Y       g DM m-2 ground
    WSOFRACTRANSLO     Amount of storage organ dry matter produced
                       by translocation of senescenced leaf organ
                       dry matter to storage organs.                   N       g DM m-2 ground
    =================  ==============================================  ======  ===========================

    **Rate variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    DLAI               Death rate of leaf area index                   Y       m2 leaf m-2 ground d-1
    DLV                Death rate of leaf dry matter                   Y       g DM m-2 d-1
    RTSUMCROPLEAFAGE   Rate of increase of crop physiological age      Y       |C|
    RWLVD              Rate of increase of dead leaf dry matter        Y       g DM m-2 d-1
    RWSOTRANSLSO       Rate of increase of of storage organ dry
                       matter
                       produced by translocation of senescenced leaf
                       organ dry matter to storage organs.             Y       g DM m-2 d-1
    =================  ==============================================  ======  ===========================

    """
    class Parameters(ParamTemplate):
        FRACTLLFENHSH = Float()
        FRACSLATB = AfgenTrait()
        FASTRANSLSO = Float()
        LAICR = Float()
        RDRB = Float()
        RDRNS = Float()
        RDRT = AfgenTrait()
        RDRSHM = Float()
        SLA_MAX = Float()
        TSUMLLIFE = Float()
        WCWET = Float()

    class RateVariables(RatesTemplate):
        DLAI = Float()
        DLV = Float()
        RTSUMCROPLEAFAGE = Float()
        RWLVD = Float()
        RWSOFASTRANSLSO = Float()
        SLA = Float()

    class StateVariables(StatesTemplate):
        TSUMCROPLEAFAGE = Float()
        WLVD = Float()
        WSOFASTRANSLSO = Float()

    def initialize(self, day, kiosk, parameters):
        TSUMCROPLEAFAGE = 0.
        WLVD = 0.
        WSOFASTRANSLSO = 0.
        self.kiosk = kiosk
        self.params = self.Parameters(parameters)
        self.rates = self.RateVariables(kiosk,
                                        publish = [
                                            "DLAI",
                                            "DLV",
                                            "RTSUMCROPLEAFAGE",
                                            "RWLVD",
                                            "RWSOFASTRANSLSO",
                                            "SLA"])
        self.states = self.StateVariables(
            kiosk,
            publish=["WLVD"],
            TSUMCROPLEAFAGE = TSUMCROPLEAFAGE,
            WLVD = WLVD,
            WSOFASTRANSLSO = WSOFASTRANSLSO
        )

    def calc_rates(self,  day, drv, delt=1):
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        # -------- AGE
        # The calculation of the physiological leaf age.
        RTSUMCROPLEAFAGE = k.DTEFF * k.EMERG - (s.TSUMCROPLEAFAGE / delt) * k.PUSHREDIST  # Deg. C

        # Relative death rate due to aging depending on leaf age and the daily average temperature.
        if s.TSUMCROPLEAFAGE - p.TSUMLLIFE >= 0:
            RDRDV = p.RDRT(drv.TEMP)
        else:
            RDRDV = 0

        # -------- SHEDDING
        # Relative death rate due to self shading, depending on a critical leaf area index at which leaf shedding is
        # induced. Leaf shedding is limited to a maximum leaf shedding per day.
        RDRSH1 = p.RDRSHM * (k.LAI-p.LAICR) / p.LAICR  # d-1

        if (RDRSH1 < 0):
            RDRSH = 0  # d-1
        elif RDRSH1 >= p.RDRSHM:
            RDRSH = p.RDRSHM  # d-1
        else:
            RDRSH = RDRSH1

        # -------- DROUGHT
        # ENSHED triggers enhanced leaf senescence due to severe drought or excessive soil water. It assumes that drought or
        # excessive water does not affect young leaves. It only affects leaves that have a reached a given fraction of the leaf
        # age.
        if k.SM - k.WCSD >= 0:
            ENHSHED1 = 0
        else:
            ENHSHED1 = 1

        if k.SM - p.WCWET >= 0:
            ENHSHED2 = 1
        else:
            ENHSHED2 = 0

        if (s.TSUMCROPLEAFAGE - p.FRACTLLFENHSH * p.TSUMLLIFE) >= 0:
            ENHSHED3 = 1
        else:
            ENHSHED3 = 0

        ENHSHED = max(ENHSHED1, ENHSHED2) * ENHSHED3

        # Relative death rate due to severe drought
        RDRSD = p.RDRB * ENHSHED  # d-1
        # --------

        # -------- NUTRIENT LIMITATION
        # Leaf death due to nutrient limitation is added on top op the relative death rate due to age, shade
        # and drought.
        RDRNS = p.RDRNS * (1-k.NPKI)  # d-1
        # --------

        # Effective relative death rate and the resulting decrease in LAI. Leaf death can only occur, when
        # the leaves are old enough.
        if s.TSUMCROPLEAFAGE - p.TSUMLLIFE >= 0:
            RDR = (max(RDRDV, RDRSH, RDRSD) + RDRNS)
        else:
            RDR = 0.

        DLAI = k.LAI * RDR * (1 - k.DORMANCY)  # m2 m-2 d-1

        # Fraction of the maximum specific leaf area index depending on the temperature sum of the crop. And its specific leaf
        # area index.
        FRACSLACROPAGE = p.FRACSLATB(k.TSUMCROP)
        SLA = p.SLA_MAX * FRACSLACROPAGE  # m2 g-1 DM

        # The rate of storage root DM production with DM supplied by the leaves before abscission.
        RWSOFASTRANSLSO = k.WLVG * RDR * p.FASTRANSLSO * (1 - k.DORMANCY)  # g storage root DM m-2 d-1

        # Decrease in leaf weight due to leaf senesence.
        DLV = k.WLVG * RDR * (1 - k.DORMANCY)  # g leaves DM m-2 d-1
        RWLVD = (DLV - RWSOFASTRANSLSO)  # g leaves DM m-2 d-1

        r.DLAI = DLAI
        r.DLV = DLV
        r.RTSUMCROPLEAFAGE = RTSUMCROPLEAFAGE
        r.RWLVD = RWLVD
        r.RWSOFASTRANSLSO = RWSOFASTRANSLSO
        r.SLA = SLA

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states
        s.WLVD += r.RWLVD
        s.WSOFASTRANSLSO += r.RWSOFASTRANSLSO
        s.TSUMCROPLEAFAGE += r.RTSUMCROPLEAFAGE