# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float

class crop_nutrient_dynamics(SimulationObject):
    """
    Class to simulate the dynamics of N, P, and K in the crop

    This class calculates the daily rates of N, P, and K uptake in the crop and their partitioning over the different
    organs.

    **Simulation parameters**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Type    Unit
    =================  ==============================================  ======  ===========================
    FR_MAX             Fraction of optimal and maximum concentration
                       of N, P, and K in each organ                    SCr     g nutrient g-1 nutrient
    K_WATER            Transpiration reduction factor at which the
                       nutrient uptake rate is reduced by a half due
                       to water stress                                 SCr     -
    NFLVD              Nitrogen concentration of dead leaves           SCr     g N g-1 DM
    KFLVD              Potassium concentration of dead leaves          SCr     g K g-1 DM
    PFLVD              Phosphorus concentration of dead leaves         SCr     g P g-1 DM
    SLOPE_NEQ_
    SOILSUPPLY_
    NEQ_PLANTUPTAKE                                                    SCr     -
    TCNPKT             Time coefficient of N, P, and K translocation   SCr     d
    =================  ==============================================  ======  ===========================

    **State variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    ANLVD              Amount of nitrogen in dead leaves               N       g N m-2 ground
    ANLVG              Amount of nitrogen in green leaves              Y       g N m-2 ground
    ANRT               Amount of nitrogen in roots                     Y       g N m-2 ground
    ANSO               Amount of nitrogen in storage organs            Y       g N m-2 ground
    ANST               Amount of nitrogen in stems                     Y       g N m-2 ground
    AKLVD              Amount of potassium in dead leaves              N       g K m-2 ground
    AKLVG              Amount of potassium in green leaves             Y       g K m-2 ground
    AKRT               Amount of potassium in roots                    Y       g K m-2 ground
    AKSO               Amount of potassium in storage organs           Y       g K m-2 ground
    AKST               Amount of potassium in stems                    Y       g K m-2 ground
    APLVD              Amount of phosphorus in dead leaves             N       g P m-2 ground
    APLVG              Amount of phosphorus in green leaves            Y       g P m-2 ground
    APRT               Amount of phosphorus in roots                   Y       g P m-2 ground
    APSO               Amount of phosphorus in storage organs          Y       g P m-2 ground
    APST               Amount of phosphorus in stems                   Y       g P m-2 ground
    =================  ==============================================  ======  ===========================

    **Rate variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    RANLVD             Rate of change of the amount of nitrogen in
                       dead leaves                                     N       g N m-2 ground d-1
    RANLVG             Rate of change of the amount of nitrogen in
                       green leaves                                    N       g N m-2 ground d-1
    RANRT              Rate of change of the amount of nitrogen in
                       roots                                           N       g N m-2 ground d-1
    RANSO              Rate of change of the amount of nitrogen in
                       storage organs                                  N       g N m-2 ground d-1
    RANST              Rate of change of the amount of nitrogen in
                       stems                                           N       g N m-2 ground d-1
    RAKLVD             Rate of change of the amount of potassium in
                       dead leaves                                     N       g K m-2 ground d-1
    RAKLVG             Rate of change of the amount of potassium in
                       green leaves                                    N       g K m-2 ground d-1
    RAKRT              Rate of change of the amount of potassium in
                       roots                                           N       g K m-2 ground d-1
    RAKSO              Rate of change of the amount of potassium in
                       storage organs                                  N       g K m-2 ground d-1
    RAKST              Rate of change of the amount of potassium in
                       stems                                           N       g K m-2 ground d-1
    RAPLVD             Rate of change of the amount of phosphorus in
                       dead leaves                                     N       g P m-2 ground d-1
    RAPLVG             Rate of change of the amount of phosphorus in
                       green leaves                                    N       g P m-2 ground d-1
    RAPRT              Rate of change of the amount of phosphorus in
                       roots                                           N       g P m-2 ground d-1
    RAPSO              Rate of change of the amount of phosphorus in
                       storage organs                                  N       g P m-2 ground d-1
    RAPST              Rate of change of the amount of phosphorus in
                       stems                                           N       g P m-2 ground d-1
    =================  ==============================================  ======  ===========================

    **Auxillary variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    RNUPTR             Daily total nitrogen uptake rate                Y       g N m-2 ground d-1
    RKUPTR             Daily total potassium uptake rate               Y       g K m-2 ground d-1
    RPUPTR             Daily total phosphorus uptake rate              Y       g P m-2 ground d-1
    =================  ==============================================  ======  ===========================

    This class is a Python implementation of the calculations related to nutrient uptake and their partitioning
    in the function nutrientdyn in the R version of the model LINTUL Cassava NPK (Adiele et al., 2022; Ezui et
    al., 2018).

    The original R version also contained calculations related to the soil supply of N, P, and K. For the
    sake of compatibility with other PCSE modules, these calculations are moved to the class
    `pcse.soil.lintul_cassava.lintul_cassava_soil_nutrient_dynamics`.

    """

    class Parameters(ParamTemplate):
        FR_MAX = Float()
        K_WATER = Float()
        SLOPE_NEQ_SOILSUPPLY_NEQ_PLANTUPTAKE = Float()
        NFLVD = Float()
        PFLVD = Float()
        KFLVD = Float()
        TCNPKT = Float()

    class RateVariables(RatesTemplate):
        RANLVG = Float()
        RANLVD = Float()
        RANST = Float()
        RANRT = Float()
        RANSO = Float()
        RAPLVG = Float()
        RAPLVD = Float()
        RAPST = Float()
        RAPRT = Float()
        RAPSO = Float()
        RAKLVG = Float()
        RAKLVD = Float()
        RAKST = Float()
        RAKRT = Float()
        RAKSO = Float()
        RNUPTR = Float()
        RPUPTR = Float()
        RKUPTR = Float()

    class StateVariables(StatesTemplate):
        ANLVG = Float()
        ANLVD = Float()
        ANST = Float()
        ANRT = Float()
        ANSO = Float()
        APLVG = Float()
        APLVD = Float()
        APST = Float()
        APRT = Float()
        APSO = Float()
        AKLVG = Float()
        AKLVD = Float()
        AKST = Float()
        AKRT = Float()
        AKSO = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        ANLVG = 0.
        ANLVD = 0.
        ANST = 0.
        ANRT = 0.
        ANSO = 0.
        APLVG = 0.
        APLVD = 0.
        APST = 0.
        APRT = 0.
        APSO = 0.
        AKLVG = 0.
        AKLVD = 0.
        AKST = 0.
        AKRT = 0.
        AKSO = 0.
        self.rates = self.RateVariables(kiosk, publish = ["RNUPTR", "RPUPTR", "RKUPTR"])
        self.states = self.StateVariables(kiosk,
                                          publish = ["ANLVG", "ANRT", "ANST", "ANSO",
                                                     "APLVG", "APRT", "APST", "APSO",
                                                     "AKLVG", "AKRT", "AKST", "AKSO"],
                                          ANLVG = ANLVG,
                                          ANLVD = ANLVD,
                                          ANST = ANST,
                                          ANRT = ANRT,
                                          ANSO = ANSO,
                                          APLVG = APLVG,
                                          APLVD = APLVD,
                                          APST = APST,
                                          APRT = APRT,
                                          APSO = APSO,
                                          AKLVG = AKLVG,
                                          AKLVD = AKLVD,
                                          AKST = AKST,
                                          AKRT = AKRT,
                                          AKSO = AKSO
                                          )

    def calc_rates(self,  day, drv, delt=1):
        # # Nutrient amounts in the crop, and the nutrient amount available for crop uptake are calculated here
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        # ---------------- Translocatable nutrient amounts
        # The amount of translocatable nutrients to the storage organs is the actual nutrient amount minus the
        # optimal amount in the plant leaves, stems and roots. For the roots it is the minimum between the
        # amount of nutrients available and as a fraction of the amount of translocatable nutrients from the
        # stem and leaves. The total translocatable nutrients is the sum of this.
        ATNLV = max(0, s.ANLVG - k.WLVG * (k.NMINLV + p.FR_MAX * (k.NMAXLV - k.NMINLV)))  # g N m-2
        ATNST = max(0, s.ANST - k.WST * (k.NMINST + p.FR_MAX * (k.NMAXST - k.NMINST)))  # g N m-2
        ATNSO = max(0, s.ANSO - k.WSO * (k.NMINSO + p.FR_MAX * (k.NMAXSO - k.NMINSO)))  # g N m-2
        ATNRT = max(0, s.ANRT - k.WRT * (k.NMINRT + p.FR_MAX * (k.NMAXRT - k.NMINRT)))  # g N m-2

        ATN = ATNLV + ATNST + ATNSO + ATNRT  # g N m-2

        ATPLV = max(0, s.APLVG - k.WLVG * (k.PMINLV + p.FR_MAX * (k.PMAXLV - k.PMINLV)))  # g P m-2
        ATPST = max(0, s.APST - k.WST * (k.PMINST + p.FR_MAX * (k.PMAXST - k.PMINST)))  # g P m-2
        ATPSO = max(0, s.APSO - k.WSO * (k.PMINSO + p.FR_MAX * (k.PMAXSO - k.PMINSO)))  # g P m-2
        ATPRT = max(0, s.APRT - k.WRT * (k.PMINRT + p.FR_MAX * (k.PMAXRT - k.PMINRT)))  # g P m-2
        ATP = ATPLV + ATPST + ATPSO + ATPRT  # g P m-2

        ATKLV = max(0, s.AKLVG - k.WLVG * (k.KMINLV + p.FR_MAX * (k.KMAXLV - k.KMINLV)))  # g K m-2
        ATKST = max(0, s.AKST - k.WST * (k.KMINST + p.FR_MAX * (k.KMAXST - k.KMINST)))  # g K m-2
        ATKSO = max(0, s.AKSO - k.WSO * (k.KMINSO + p.FR_MAX * (k.KMAXSO - k.KMINSO)))  # g K m-2
        ATKRT = max(0, s.AKRT - k.WRT * (k.KMINRT + p.FR_MAX * (k.KMAXRT - k.KMINRT)))  # g K m-2
        ATK = ATKLV + ATKST + ATKSO + ATKRT  # g K m-2
        # # --------------
        #
        # # ---------------- Nutrient demand
        # # The nutrient demand is calculated as the difference between the amount of translocatable nutrients and
        # # the maximum nutrient content. The total nutrient demand is the sum of the demands of the different organs.
        #
        NDEML = max(k.NMAXLV * k.WLVG - s.ANLVG, 0)  # g N m-2
        NDEMS = max(k.NMAXST * k.WST - s.ANST, 0)  # g N m-2
        NDEMR = max(k.NMAXRT * k.WRT - s.ANRT, 0)  # g N m-2
        NDEMSO = max(k.NMAXSO * k.WSO - s.ANSO, 0)  # g N m-2
        NDEMTO = max(0, (NDEML + NDEMS + NDEMSO + NDEMR))  # g N m-2

        PDEML = max(k.PMAXLV * k.WLVG - s.APLVG, 0)  # g P m-2
        PDEMS = max(k.PMAXST * k.WST - s.APST, 0)  # g P m-2
        PDEMR = max(k.PMAXRT * k.WRT - s.APRT, 0)  # g P m-2
        PDEMSO = max(k.PMAXSO * k.WSO - s.APSO, 0)  # g P m-2
        PDEMTO = max(0, (PDEML + PDEMS + PDEMSO + PDEMR))  # g P m-2

        KDEML = max(k.KMAXLV * k.WLVG - s.AKLVG, 0)  # g K m-2
        KDEMS = max(k.KMAXST * k.WST - s.AKST, 0)  # g K m-2
        KDEMR = max(k.KMAXRT * k.WRT - s.AKRT, 0)  # g K m-2
        KDEMSO = max(k.KMAXSO * k.WSO - s.AKSO, 0)  # g K m-2
        KDEMTO = max(0, (KDEML + KDEMS + KDEMSO + KDEMR))  # g K m-2
        # ---------------

        # --------------- Net nutrient translocation in the crop
        # Internal relocation of nutrients depends on relative content
        # When there is no demand, there should also be no redistribution

        # Computation of the translocation of nutrients from the different organs
        # Daily redistribution to balance nutrient contents between all organs
        # Based on relative demand of organs
        if p.TCNPKT * NDEMTO == 0:
            RNTLV = 0.
            RNTST = 0.
            RNTSO = 0.
            RNTRT = 0.
        else:
            RNTLV = ((NDEML / NDEMTO) * ATN - ATNLV) / p.TCNPKT
            RNTST = ((NDEMS / NDEMTO) * ATN - ATNST) / p.TCNPKT
            RNTSO = ((NDEMSO / NDEMTO) * ATN - ATNSO) / p.TCNPKT
            RNTRT = ((NDEMR / NDEMTO) * ATN - ATNRT) / p.TCNPKT

        if p.TCNPKT * PDEMTO == 0:
            RPTLV = 0.
            RPTST = 0.
            RPTSO = 0.
            RPTRT = 0.
        else:
            RPTLV = ((PDEML / PDEMTO) * ATP - ATPLV) / p.TCNPKT
            RPTST = ((PDEMS / PDEMTO) * ATP - ATPST) / p.TCNPKT
            RPTSO = ((PDEMSO / PDEMTO) * ATP - ATPSO) / p.TCNPKT
            RPTRT = ((PDEMR / PDEMTO) * ATP - ATPRT) / p.TCNPKT

        if p.TCNPKT * KDEMTO == 0:
            RKTLV = 0.
            RKTST = 0.
            RKTSO = 0.
            RKTRT = 0.
        else:
            RKTLV = ((KDEML / KDEMTO) * ATK - ATKLV) / p.TCNPKT
            RKTST = ((KDEMS / KDEMTO) * ATK - ATKST) / p.TCNPKT
            RKTSO = ((KDEMSO / KDEMTO) * ATK - ATKSO) / p.TCNPKT
            RKTRT = ((KDEMR / KDEMTO) * ATK - ATKRT) / p.TCNPKT

        # # ---------------
        TINY = 10E-9
        if (abs(RNTLV + RNTST + RNTSO + RNTRT) > TINY):
            print("UNRELIABLE RESULTS!! Internal N reallocation must be net 0")
            print(f"RNTLV = {RNTLV}, RNTST = {RNTST}, RNTSO = {RNTSO}, RNTRT = {RNTRT}")
        if (abs(RPTLV + RPTST + RPTSO + RPTRT) > TINY):
            print("UNRELIABLE RESULTS!! Internal P reallocation must be net 0")
            print(f"RPTLV = {RPTLV}, RPTST = {RPTST}, RPTSO = {RPTSO}, RPTRT = {RPTRT}")
        if (abs(RKTLV + RKTST + RKTSO + RKTRT) > TINY):
            print("UNRELIABLE RESULTS!! Internal K reallocation must be net 0")
            print(f"RKTLV = {RKTLV}, RKTST = {RKTST}, RKTSO = {RKTSO}, RKTRT = {RKTRT}")

        # --------------- Nutrient uptake
        # Nutrient uptake from the soil depends on the soil moisture. It is assumed
        # that nutrient uptake reduces monod-like with growth rate reduction due to:
        # 1) Low soil water supply: uptake rates are 50% when soil supply rates equals max-uptake rates
        # 2) Low soil nutrient supply: uptake rates are 50% when max soil supply rates equals max-uptake rates
        WLIMIT = k.RFTRA / (p.K_WATER + k.RFTRA)

        # Maximum amounts of nutrients for the given amount of biomass
        NMAX = k.NMAXLV * k.WLVG + k.NMAXST * k.WST + k.NMAXRT * k.WRT + k.NMAXSO * k.WSO
        PMAX = k.PMAXLV * k.WLVG + k.PMAXST * k.WST + k.PMAXRT * k.WRT + k.PMAXSO * k.WSO
        KMAX = k.KMAXLV * k.WLVG + k.KMAXST * k.WST + k.KMAXRT * k.WRT + k.KMAXSO * k.WSO

        # Nutrient equivalents in the soil and maximum uptake of equivalents based on optimum ratios
        # g m-2
        if PMAX == 0:
            NMAX2PMAX = 0
        else:
            NMAX2PMAX = NMAX / PMAX
        if KMAX == 0:
            NMAX2KMAX = 0
        else:
            NMAX2KMAX = NMAX / KMAX

        NUTEQ_SOIL = k.NAVAIL + k.PAVAIL * NMAX2PMAX + k.KAVAIL * NMAX2KMAX

        # g m-2
        NUTEQ_DEMAND = NDEMTO + PDEMTO * NMAX2PMAX + KDEMTO * NMAX2KMAX

        # The parameter SLOPE_NEQ_SOIL_PEQUPTAKE determines the ratio between biomass and nutrients
        # Plant uptake is proportional to nutrient supply under optimum ratios.
        # uptake rates increase when soil supply is larger.
        # Interaction effects are not accounted for here.

        # If uptake isn't adequate, concentrations will decrease first, later growth rates decrease.
        # gNEQ m-2 d-1     g m-2           d-1,       d-1                              g m-2
        RMAX_UPRATE = min(NUTEQ_DEMAND / delt, p.SLOPE_NEQ_SOILSUPPLY_NEQ_PLANTUPTAKE * NUTEQ_SOIL)

        # If actual ratios are suboptimal, uptake of excess nutrients is relatively increased.
        # A suboptimal ratio in the soil results in a suboptimal ratio in the plant when supply is smaller than demand.
        # For example, a soil full of N gives large N uptake, but relatively small P and K uptake.
        # Actual uptake of N is limited by maximum N contensts, reducing uptake of P and K contents in the plant
        # to suboptimal levels. In extremis, oversupply of one nutrient can reduce uptake of another.
        #
        # The fraction of NEQ taken up that is determined by N equals: NEQ_M =  NMINT/NUEQ_SOIL
        # The fraction of NEQ taken up that is determined by P equals: NEQ_P = (NMAX/PMAX)*PMINT/NUTEQ_SOIL
        # The fraction of NEQ taken up that is determined by K equals: NEQ_K = (NMAX/KMAX)*KMINT/NUTEQ_SOIL
        # To translate NEQ_P uptake back to P uptake, this needs to be multiplied by:
        # PUP = NEQ_P * PMAX / NMAX
        # PUK = NEQ_K * KMAX / NMAX

        # The fraction of NEQ taken up that is determined by N equals: NMINT/NUEQ_SOIL
        if NUTEQ_SOIL == 0:
            RNUPTR = 0.
            RPUPTR = 0.
            RKUPTR = 0
        else:
            # gN m-2 d-1 = gNEQ m-2 d-1 * gN m-2 * gNEQ-1 m2
            RNUPTR = (RMAX_UPRATE * k.NAVAIL / NUTEQ_SOIL) * WLIMIT

            # gP m-2 d-1 = gNEQ m-2 d-1 gN kgDM-1 * kgDM gP-1 * gP m-2 * gNEQ-1 m2
            RPUPTR = (RMAX_UPRATE * k.PAVAIL / NUTEQ_SOIL) * WLIMIT

            # gK m-2 d-1 = gNEQ m-2 d-1 gN kgDM-1 * kgDM gK-1 * gK m-2 * gNEQ-1 m2
            RKUPTR = (RMAX_UPRATE * k.KAVAIL / NUTEQ_SOIL) * WLIMIT

        # Actual uptake is limited by demand based on maximum concentrations for standing biomass
        # Uptake rate should not exceed maximum soil supply.
        RNUPTR = min(k.NAVAIL / delt, RNUPTR, NDEMTO / delt)
        RPUPTR = min(k.PAVAIL / delt, RPUPTR, PDEMTO / delt)
        RKUPTR = min(k.KAVAIL / delt, RKUPTR, KDEMTO / delt)
        # -------------

        # ------------- Partitioning
        # to compute the partitioning of the total N/P/K uptake rates (NUPTR, PUPTR, KUPTR)
        # over the leaves, stem, and roots (kg N/P/K ha-1 d-1)
        # concentrations are balanced, so distribution based on weight proportions
        WTOT = k.WLVG + k.WST + k.WSO + k.WRT

        if WTOT == 0:
            RNULV = 0.
            RNUST = 0.
            RNUSO = 0.
            RNURT = 0.

            RPULV = 0.
            RPUST = 0.
            RPUSO = 0.
            RPURT = 0.

            RKULV = 0.
            RKUST = 0.
            RKUSO = 0.
            RKURT = 0.
        else:
            RNULV = (k.WLVG / WTOT) * RNUPTR
            RNUST = (k.WST / WTOT) * RNUPTR
            RNUSO = (k.WSO / WTOT) * RNUPTR
            RNURT = (k.WRT / WTOT) * RNUPTR

            RPULV = (k.WLVG / WTOT) * RPUPTR
            RPUST = (k.WST / WTOT) * RPUPTR
            RPUSO = (k.WSO / WTOT) * RPUPTR
            RPURT = (k.WRT / WTOT) * RPUPTR

            RKULV = (k.WLVG / WTOT) * RKUPTR
            RKUST = (k.WST / WTOT) * RKUPTR
            RKUSO = (k.WSO / WTOT) * RKUPTR
            RKURT = (k.WRT / WTOT) * RKUPTR
        #
        # ------------ Nutrient redistribution because of cutting.
        # A negative rate in e.g. RNCUTTING adds to other components.
        RANCUTLV = -k.RNCUTTING * k.FLV  # g N m-2 d-1
        RANCUTST = -k.RNCUTTING * k.FST  # g N m-2 d-1
        RANCUTRT = -k.RNCUTTING * k.FRT  # g N m-2 d-1
        RANCUTSO = -k.RNCUTTING * k.FSO  # g N m-2 d-1

        RAPCUTLV = -k.RPCUTTING * k.FLV  # g P m-2 d-1
        RAPCUTST = -k.RPCUTTING * k.FST  # g P m-2 d-1
        RAPCUTRT = -k.RPCUTTING * k.FRT  # g P m-2 d-1
        RAPCUTSO = -k.RPCUTTING * k.FSO  # g P m-2 d-1

        RAKCUTLV = -k.RKCUTTING * k.FLV  # g K m-2 d-1
        RAKCUTST = -k.RKCUTTING * k.FST  # g K m-2 d-1
        RAKCUTRT = -k.RKCUTTING * k.FRT  # g K m-2 d-1
        RAKCUTSO = -k.RKCUTTING * k.FSO  # g K m-2 d-1
        # ------------

        # ------------ Nutrient redistribution because of leaf death
        if k.WLVG > 0:
            # Nutrients lost due to dying leaves
            RANLVD = k.RWLVD * p.NFLVD  # g N m-2 d-1
            RAPLVD = k.RWLVD * p.PFLVD  # g P m-2 d-1
            RAKLVD = k.RWLVD * p.KFLVD  # g K m-2 d-1
            # Total nutrients in dying leaves
            RNDLVG = k.RWLVD * (s.ANLVG / k.WLVG)  # g N m-2 d-1
            RPDLVG = k.RWLVD * (s.APLVG / k.WLVG)  # g P m-2 d-1
            RKDLVG = k.RWLVD * (s.AKLVG / k.WLVG)  # g K m-2 d-1
        else:
            RNDLVG = 0.
            RPDLVG = 0.
            RKDLVG = 0.
            RANLVD = 0.
            RAPLVD = 0.
            RAKLVD = 0.

        # What is not lost to dead leaves most be redistributed to other organs
        RNDLV_REDIST = max(0, RNDLVG - RANLVD)  # g N m-2 d-1
        RPDLV_REDIST = max(0, RPDLVG - RAPLVD)
        RKDLV_REDIST = max(0, RKDLVG - RAKLVD)

        # ------------

        # ------------ Nutrient redistribution because of storage root DM redistribution after dormancy
        # DM to the leaves, with new at maximum NPK concentrations
        #             g DM m-2 d-1 * (gN m-2 d-1 * gDM-1 m2 d)
        RANSO2LVLV = k.RREDISTLVG * k.NMAXLV * k.PUSHREDIST  # g N m-2 d-1
        RAPSO2LVLV = k.RREDISTLVG * k.PMAXLV * k.PUSHREDIST  # g P m-2 d-1
        RAKSO2LVLV = k.RREDISTLVG * k.KMAXLV * k.PUSHREDIST  # g K m-2 d-1

        # DM loss of the storage roots
        if k.WSO == 0:
            RANSO2LVSO = 0.
            RAPSO2LVSO = 0.
            RAKSO2LVSO = 0.
        else:
            RANSO2LVSO = k.RREDISTSO * (s.ANSO / k.WSO)
            RAPSO2LVSO = k.RREDISTSO * (s.APSO / k.WSO)
            RAKSO2LVSO = k.RREDISTSO * (s.AKSO / k.WSO)

        # ------------- Rate of change of N/P/K in crop organs
        #        uptake + net translocation + cutting
        # N relocated to stem, P+K to storate roots
        RANLVG = RNULV + RNTLV + RANCUTLV + RANSO2LVLV - RNDLVG  # g N m-2 d-1
        RANST = RNUST + RNTST + RANCUTST + RNDLV_REDIST  # g N m-2 d-1
        RANRT = RNURT + RNTRT + RANCUTRT  # g N m-2 d-1
        RANSO = RNUSO + RNTSO + RANCUTSO - RANSO2LVSO  # g N m-2 d-1

        RAPLVG = RPULV + RPTLV + RAPCUTLV + RAPSO2LVLV - RPDLVG  # g P m-2 d-1
        RAPST = RPUST + RPTST + RAPCUTST  # g P m-2 d-1
        RAPRT = RPURT + RPTRT + RAPCUTRT  # g P m-2 d-1
        RAPSO = RPUSO + RPTSO + RAPCUTSO - RAPSO2LVSO + RPDLV_REDIST  # g P m-2 d-1

        RAKLVG = RKULV + RKTLV + RAKCUTLV + RAKSO2LVLV - RKDLVG  # g K m-2 d-1
        RAKST = RKUST + RKTST + RAKCUTST  # g K m-2 d-1
        RAKRT = RKURT + RKTRT + RAKCUTRT  # g K m-2 d-1
        RAKSO = RKUSO + RKTSO + RAKCUTSO - RAKSO2LVSO + RKDLV_REDIST  # g K m-2 d-1  #

        if (k.EMERG == 1) & (k.WST == 0):
            RANLVG = 0.
            RANLVD = 0.
            RANST = 0.
            RANRT = 0.
            RANSO = 0.

            RAPLVG = 0.
            RAPLVD = 0.
            RAPST = 0.
            RAPRT = 0.
            RAPSO = 0.

            RAKLVG = 0.
            RAKLVD = 0.
            RAKST = 0.
            RAKRT = 0.
            RAKSO = 0.
        else:
            pass

        r.RANLVG = RANLVG
        r.RANLVD = RANLVD
        r.RANST = RANST
        r.RANRT = RANRT
        r.RANSO = RANSO
        r.RAPLVG = RAPLVG
        r.RAPLVD = RAPLVD
        r.RAPST = RAPST
        r.RAPRT = RAPRT
        r.RAPSO = RAPSO
        r.RAKLVG = RAKLVG
        r.RAKLVD = RAKLVD
        r.RAKST = RAKST
        r.RAKRT = RAKRT
        r.RAKSO = RAKSO
        r.RNUPTR = RNUPTR
        r.RPUPTR = RPUPTR
        r.RKUPTR = RKUPTR

    def integrate(self, day, drv, delt = 1):
        r = self.rates
        s = self.states

        s.ANLVG += delt * r.RANLVG
        s.ANLVD += delt * r.RANLVD
        s.ANST += delt * r.RANST
        s.ANRT += delt * r.RANRT
        s.ANSO += delt * r.RANSO
        s.APLVG += delt * r.RAPLVG
        s.APLVD += delt * r.RAPLVD
        s.APST += delt * r.RAPST
        s.APRT += delt * r.RAPRT
        s.APSO += delt * r.RAPSO
        s.AKLVG += delt * r.RAKLVG
        s.AKLVD += delt * r.RAKLVD
        s.AKST += delt * r.RAKST
        s.AKRT += delt * r.RAKRT
        s.AKSO += delt * r.RAKSO