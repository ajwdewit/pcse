# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

import numpy as np

from pcse.traitlets import Bool, Float, Instance
from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.crop.lintul_cassava.lintul_cassava_util import afgen2cols

class npk_stress(SimulationObject):
    """
    Class to simulate nutrient stress indices

    This class calculates nutrient stress indices that are used to calculate how deficiencies of N, P, and K affect
    the dry matter production, leaf senescence and the the biomass partitioning.

    **Simulation parameters**

    =================  ================================================  ======  ===========================
    Name               Description                                       Type    Unit
    =================  ================================================  ======  ===========================
    FR_MAX             Fraction of optimal and maximum concentration
                       of nitrogen, potassimu, and phosphorus in
                       each organ                                         SCr     g nutrient g-1 nutrient
    NMINMAXLV          Minimum and maximum N concentrations in leaves
                       as a function of temperature sum.                  TCr     g N kg-1 DM
    NMINMAXRT          Minimum and maximum N concentrations in roots
                       as a function of temperature sum.                  TCr     g N kg-1 DM
    NMINMAXSO          Minimum and maximum N concentrations in storage
                       organs as a function of temperature sum.           TCr     g N kg-1 DM
    NMINMAXST          Minimum and maximum N concentrations in stems
                       as a function of temperature sum.                  TCr     g N kg-1 DM
    KMINMAXLV          Minimum and maximum P concentrations in leaves
                       as a function of temperature sum.                  TCr     g K kg-1 DM
    K_MAX              Maximum value of K in Monod function.              SCr     -
    K_NPK_NI           K value in Monod relationship to reduce
                       influence of slightly lower NI value. A higher
                       values give quicker stress.                        SCr     -
    KMINMAXRT          Minimum and maximum K concentrations in roots
                       as a function of temperature sum.                  TCr     g K kg-1 DM
    KMINMAXSO          Minimum and maximum K concentrations in storage
                       organs as a function of temperature sum.           TCr     g K kg-1 DM
    KMINMAXST          Minimum and maximum K concentrations in stems
                       as a function of temperature sum.                  TCr     g K kg-1 DM
    PMINMAXLV          Minimum and maximum K concentrations in leaves
                       as a function of temperature sum.                  TCr     g P kg-1 DM
    PMINMAXRT          Minimum and maximum P concentrations in roots
                       as a function of temperature sum.                  TCr     g P kg-1 DM
    PMINMAXSO          Minimum and maximum P concentrations in storage
                       organs as a function of temperature sum.           TCr     g P kg-1 DM
    PMINMAXST          Minimum and maximum P concentrations in stems
                       as a function of temperature sum.                  TCr     g P kg-1 DM
    TSUM_NPKI          Temperature sum below which there is
                       no nutrient stress                                 SCr     |C| d
    =================  ================================================  ======  ===========================

    **Auxillary variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    NNI                Nitrogen Nutrition Index                        Y       -
    PNI                Phosphorus Nutrition Index                      Y       -
    KNI                Potassium Nutrition Index                       Y       -
    NPKI               Lumped N,P, and K Nutrition Index               Y       -

    NMINLV             Minimum N concentration in leaves               Y       g N kg-1 DM
    PMINLV             Minimum P concentration in leaves               Y       g P kg-1 DM
    KMINLV             Minimum K concentration in leaves               Y       g K kg-1 DM
    NMINST             Minimum N concentration in stems                Y       g N kg-1 DM
    PMINST             Minimum P concentration in stems                Y       g P kg-1 DM
    KMINST             Minimum K concentration in stems                Y       g K kg-1 DM
    NMINSO             Minimum N concentration in storage organs       Y       g N kg-1 DM
    PMINSO             Minimum P concentration in storage organs       Y       g P kg-1 DM
    KMINSO             Minimum K concentration in storage organs       Y       g K kg-1 DM
    NMINRT             Minimum N concentration in roots                Y       g N kg-1 DM
    PMINRT             Minimum K concentration in roots                Y       g N kg-1 DM
    KMINRT             Minimum P concentration in roots                Y       g N kg-1 DM
    NMAXLV             Maximum N concentration in leaves               Y       g N kg-1 DM
    PMAXLV             Maximum P concentration in leaves               Y       g P kg-1 DM
    KMAXLV             Maximum K concentration in leaves               Y       g K kg-1 DM
    NMAXST             Maximum N concentration in stems                Y       g N kg-1 DM
    PMAXST             Maximum P concentration in stems                Y       g P kg-1 DM
    KMAXST             Maximum K concentration in stems                Y       g K kg-1 DM
    NMAXSO             Maximum N concentration in storage organs       Y       g N kg-1 DM
    PMAXSO             Maximum P concentration in storage organs       Y       g P kg-1 DM
    KMAXSO             Maximum K concentration in storage organs       Y       g K kg-1 DM
    NMAXRT             Maximum N concentration in roots                Y       g N kg-1 DM
    PMAXRT             Maximum K concentration in roots                Y       g N kg-1 DM
    KMAXRT             Maximum P concentration in roots                Y       g N kg-1 DM
    =================  ==============================================  ======  ===========================
    """

    NUTRIENT_LIMITED = True

    class Parameters(ParamTemplate):
        FR_MAX = Float()
        K_MAX = Float()
        K_NPK_NI = Float()
        TSUM_NPKI = Float()
        NMINMAXLV = Instance(list)
        PMINMAXLV = Instance(list)
        KMINMAXLV = Instance(list)
        NMINMAXRT = Instance(list)
        PMINMAXRT = Instance(list)
        KMINMAXRT = Instance(list)
        NMINMAXST = Instance(list)
        PMINMAXST = Instance(list)
        KMINMAXST = Instance(list)
        NMINMAXSO = Instance(list)
        PMINMAXSO = Instance(list)
        KMINMAXSO = Instance(list)

    class RateVariables(RatesTemplate):
        NNI = Float()
        PNI = Float()
        KNI = Float()
        NPKI = Float()

        NMINLV = Float()
        PMINLV = Float()
        KMINLV = Float()
        NMINST = Float()
        PMINST = Float()
        KMINST = Float()
        NMINSO = Float()
        PMINSO = Float()
        KMINSO = Float()
        NMINRT = Float()
        PMINRT = Float()
        KMINRT = Float()
        NMAXLV = Float()
        PMAXLV = Float()
        KMAXLV = Float()
        NMAXST = Float()
        PMAXST = Float()
        KMAXST = Float()
        NMAXSO = Float()
        PMAXSO = Float()
        KMAXSO = Float()
        NMAXRT = Float()
        PMAXRT = Float()
        KMAXRT = Float()

    class StateVariables(StatesTemplate):
        pass

    def initialize(self, day, kiosk, parameters):
        self.kiosk = kiosk
        self.params = self.Parameters(parameters)
        self.rates = self.RateVariables(kiosk,
                                        publish = ["NPKI", "NNI", "PNI", "KNI",
                                                   "NMINLV", "PMINLV","KMINLV",
                                                   "NMINST", "PMINST","KMINST",
                                                   "NMINSO", "PMINSO","KMINSO",
                                                   "NMINRT", "PMINRT", "KMINRT",
                                                   "NMAXLV", "PMAXLV", "KMAXLV",
                                                   "NMAXST", "PMAXST", "KMAXST",
                                                   "NMAXSO", "PMAXSO", "KMAXSO",
                                                   "NMAXRT", "PMAXRT", "KMAXRT"])
        self.states = self.StateVariables(
            kiosk,
            publish=[]
        )

    def __call__(self, day, drv, delt = 1):
        k = self.kiosk
        p = self.params
        r = self.rates

        # The nutrient limitation is based on the nutrient concentrations in the organs of the crop. A nutrition index
        # is calculated to quantify nutrient limitation.

        # Minimum and maximum nutrient concentrations in the leaves
        NMINLV = afgen2cols(p.NMINMAXLV, k.TSUMCROP, 1)  # g N g-1 DM
        PMINLV = afgen2cols(p.PMINMAXLV, k.TSUMCROP, 1)  # g N g-1 DM
        KMINLV = afgen2cols(p.KMINMAXLV, k.TSUMCROP, 1)  # g N g-1 DM
        NMAXLV = afgen2cols(p.NMINMAXLV, k.TSUMCROP, 2)  # g N g-1 DM
        PMAXLV = afgen2cols(p.PMINMAXLV, k.TSUMCROP, 2)  # g N g-1 DM
        KMAXLV = afgen2cols(p.KMINMAXLV, k.TSUMCROP, 2)  # g N g-1 DM

        # # Minimum and maximum concentrations in the stems
        NMINST = afgen2cols(p.NMINMAXST, k.TSUMCROP, 1)  # g N g-1 DM
        PMINST = afgen2cols(p.PMINMAXST, k.TSUMCROP, 1)  # g N g-1 DM
        KMINST = afgen2cols(p.KMINMAXST, k.TSUMCROP, 1)  # g N g-1 DM
        NMAXST = afgen2cols(p.NMINMAXST, k.TSUMCROP, 2)  # g N g-1 DM
        PMAXST = afgen2cols(p.PMINMAXST, k.TSUMCROP, 2)  # g N g-1 DM
        KMAXST = afgen2cols(p.KMINMAXST, k.TSUMCROP, 2)  # g N g-1 DM

        # # Minimum and maximum nutrient concentrations in the storage organs
        NMINSO = afgen2cols(p.NMINMAXSO, k.TSUMCROP, 1)  # g N g-1 DM
        PMINSO = afgen2cols(p.PMINMAXSO, k.TSUMCROP, 1)  # g N g-1 DM
        KMINSO = afgen2cols(p.KMINMAXSO, k.TSUMCROP, 1)  # g N g-1 DM
        NMAXSO = afgen2cols(p.NMINMAXSO, k.TSUMCROP, 2)  # g N g-1 DM
        PMAXSO = afgen2cols(p.PMINMAXSO, k.TSUMCROP, 2)  # g N g-1 DM
        KMAXSO = afgen2cols(p.KMINMAXSO, k.TSUMCROP, 2)  # g N g-1 DM

        # # Minimum and maximum nutrient concentrations in the roots
        NMINRT = afgen2cols(p.NMINMAXRT, k.TSUMCROP, 1)  # g N g-1 DM
        PMINRT = afgen2cols(p.PMINMAXRT, k.TSUMCROP, 1)  # g N g-1 DM
        KMINRT = afgen2cols(p.KMINMAXRT, k.TSUMCROP, 1)  # g N g-1 DM
        NMAXRT = afgen2cols(p.NMINMAXRT, k.TSUMCROP, 2)  # g N g-1 DM
        PMAXRT = afgen2cols(p.PMINMAXRT, k.TSUMCROP, 2)  # g N g-1 DM
        KMAXRT = afgen2cols(p.KMINMAXRT, k.TSUMCROP, 2)  # g N g-1 DM

        # ---------------- Nutrient concentrations
        # Minimum nutrient content in the living biomass
        NMIN = k.WLVG * NMINLV + k.WST * NMINST + k.WSO * NMINSO  # g N m-2
        PMIN = k.WLVG * PMINLV + k.WST * PMINST + k.WSO * PMINSO  # g P m-2
        KMIN = k.WLVG * KMINLV + k.WST * KMINST + k.WSO * KMINSO  # g K m-2

        # Maximum nutrient content in the living biomass
        NMAX = NMAXLV * k.WLVG + NMAXST * k.WST + NMAXSO * k.WSO  # g N m-2
        PMAX = PMAXLV * k.WLVG + PMAXST * k.WST + PMAXSO * k.WSO  # g P m-2
        KMAX = KMAXLV * k.WLVG + KMAXST * k.WST + KMAXSO * k.WSO  # g K m-2

        # Optimal nutrient content in the living biomass
        NOPT = NMIN + p.FR_MAX * (NMAX - NMIN)  # g N m-2
        POPT = PMIN + p.FR_MAX * (PMAX - PMIN)  # g P m-2
        KOPT = KMIN + p.FR_MAX * (KMAX - KMIN)  # g K m-2
        # ----------------

        # ---------------- Actual nutrient concentrations

        # Actual nutrient amounts in the living biomass
        NACT = k.ANLVG + k.ANST + k.ANSO  # g N m-2
        PACT = k.APLVG + k.APST + k.APSO  # g P m-2
        KACT = k.AKLVG + k.AKST + k.AKSO  # g K m-2

        # -------------- Nutrition Indices
        if NOPT - NMIN == 0:
            NNI = 0
        else:
            NNI = (NACT - NMIN) / (NOPT - NMIN)
        if POPT - PMIN == 0:
            PNI = 0
        else:
            PNI = (PACT - PMIN) / (POPT - PMIN)
        if KOPT - KMIN == 0:
            KNI = 0
        else:
            KNI = (KACT - KMIN) / (KOPT - KMIN)

        NNI = min(1, max(0, NNI))
        PNI = min(1, max(0, PNI))
        KNI = min(1, max(0, KNI))

        # Combined effect.
        # Multiplication allows to have extra growth reduction if multiple nutrients are deficient
        # The "Monod" acts as scalar to reduce effect of minor deficiencies that do not affect growth rates
        # but are compensated by dilution.
        # A mirrored Monod function to determine effect of N, P and K stress on NPKI
        NPKI = self.Mirrored_Monod(x=NNI * PNI * KNI, K=p.K_NPK_NI, Kmax=p.K_MAX)

        # Nutrient limitation reduction factor when nutrient limition is switched on
        if self.NUTRIENT_LIMITED:
            # Simple based on daily values
            NPKI = max(0, min(1, NPKI))  # (-)
            # Shortly after emergence nutrient stress does not occur
            if k.TSUMCROP < p.TSUM_NPKI:
                NPKI = 1
            else:
                pass
        else:
            NPKI =1

        r.NNI = NNI
        r.PNI = PNI
        r.KNI = KNI
        r.NPKI = NPKI

        r.NMINLV = NMINLV
        r.PMINLV = PMINLV
        r.KMINLV = KMINLV
        r.NMINST = NMINST
        r.PMINST = PMINST
        r.KMINST = KMINST
        r.NMINSO = NMINSO
        r.PMINSO = PMINSO
        r.KMINSO = KMINSO
        r.NMINRT = NMINRT
        r.PMINRT = PMINRT
        r.KMINRT = KMINRT
        r.NMAXLV = NMAXLV
        r.PMAXLV = PMAXLV
        r.KMAXLV = KMAXLV
        r.NMAXST = NMAXST
        r.PMAXST = PMAXST
        r.KMAXST = KMAXST
        r.NMAXSO = NMAXSO
        r.PMAXSO = PMAXSO
        r.KMAXSO = KMAXSO
        r.NMAXRT = NMAXRT
        r.PMAXRT = PMAXRT
        r.KMAXRT = KMAXRT

    def Mirrored_Monod(self, x, K, Kmax=4):
        if K <= Kmax:
            C = K + 1
            if x == 0:
                y = 0
            else:
                y = C * x / (x + K)
        else:
            K = max(0, 2 * Kmax - K)
            C = K + 1
            x = 1 - x
            if x == 0:
                y = 0
            else:
                y = C * x / (x + K)
            y = 1 - y
        return y