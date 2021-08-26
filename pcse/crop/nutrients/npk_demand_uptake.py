# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
from collections import namedtuple

from ...base import StatesTemplate, ParamTemplate, SimulationObject, RatesTemplate
from ...decorators import prepare_rates, prepare_states
from ...traitlets import HasTraits, Float, Int, Instance
from ...util import AfgenTrait

MaxNutrientConcentrations = namedtuple("MaxNutrientConcentrations",
                                       ["NMAXLV", "PMAXLV", "KMAXLV",
                                        "NMAXST", "PMAXST", "KMAXST",
                                        "NMAXRT", "PMAXRT", "KMAXRT",
                                        "NMAXSO", "PMAXSO", "KMAXSO"])

class NPK_Demand_Uptake(SimulationObject):
    """Calculates the crop N/P/K demand and its uptake from the soil.

    Crop N/P/K demand is calculated as the difference between the
    actual N/P/K concentration (kg N/P/K per kg biomass) in the
    vegetative plant organs (leaves, stems and roots) and the maximum
    N/P/K concentration for each organ. N/P/K uptake is then estimated
    as the minimum of supply from the soil and demand from the crop.

    Nitrogen fixation (leguminous plants) is calculated by assuming that a
    fixed fraction of the daily N demand is supplied by nitrogen fixation.
    The remaining part has to be supplied by the soil.

    The N/P/K demand of the storage organs is calculated in a somewhat
    different way because it is assumed that the demand from the storage
    organs is fulfilled by translocation of N/P/K from the leaves, stems
    and roots. So Therefore the uptake of the storage organs is calculated
    as the minimum of the translocatable N/P/K (supply) and the demand from
    the storage organs. Moreover, there is time coefficient for translocation
    which takes into account that there is a delay in the availability of
    translocatable N/P/K

    **Simulation parameters**

    ============  =============================================  ======================
     Name          Description                                    Unit
    ============  =============================================  ======================
    NMAXLV_TB      Maximum N concentration in leaves as          kg N kg-1 dry biomass
                   function of DVS
    PMAXLV_TB      As for P                                      kg P kg-1 dry biomass
    KMAXLV_TB      As for K                                      kg K kg-1 dry biomass

    NMAXRT_FR      Maximum N concentration in roots as fraction  -
                   of maximum N concentration in leaves
    PMAXRT_FR      As for P                                      -
    KMAXRT_FR      As for K                                      -

    NMAXST_FR      Maximum N concentration in stems as fraction  -
                   of maximum N concentration in leaves
    PMAXST_FR      As for P                                      -
    KMAXST_FR      As for K                                      -

    NMAXSO         Maximum N concentration in storage organs     kg N kg-1 dry biomass
    PMAXSO         As for P                                      kg P kg-1 dry biomass
    KMAXSO         As for K                                      kg K kg-1 dry biomass

    NCRIT_FR       Critical N concentration as fraction of       -
                   maximum N concentration for vegetative
                   plant organs as a whole (leaves + stems)
    PCRIT_FR       As for P                                      -
    KCRIT_FR       As for K                                      -

    TCNT           Time coefficient for N translation to         days
                   storage organs
    TCPT           As for P                                      days
    TCKT           As for K                                      days

    NFIX_FR        fraction of crop nitrogen uptake by           kg N kg-1 dry biomass
                   biological fixation
    DVS_NPK_STOP   Development stage after which no nutrients    -
                   are taken up from the soil by the crop.
    ============  =============================================  ======================

    **State variables**

    ==========  ================================================= ==== ============
     Name        Description                                      Pbl      Unit
    ==========  ================================================= ==== ============
    NdemandLV     N Demand in living leaves                         N   |kg N ha-1|
    NdemandST     N Demand in living stems                          N   |kg N ha-1|
    NdemandRT     N Demand in living roots                          N   |kg N ha-1|
    NdemandSO     N Demand in storage organs                        N   |kg N ha-1|

    PdemandLV     P Demand in living leaves                         N   |kg P ha-1|
    PdemandST     P Demand in living stems                          N   |kg P ha-1|
    PdemandRT     P Demand in living roots                          N   |kg P ha-1|
    PdemandSO     P Demand in storage organs                        N   |kg P ha-1|

    KdemandLV     K Demand in living leaves                         N   |kg K ha-1|
    KdemandST     K Demand in living stems                          N   |kg K ha-1|
    KdemandRT     K Demand in living roots                          N   |kg K ha-1|
    KdemandSO     K Demand in storage organs                        N   |kg K ha-1|
    ==========  ================================================= ==== ============


    **Rate variables**

    ===========  ================================================= ==== ================
     Name         Description                                      Pbl      Unit
    ===========  ================================================= ==== ================
    RNuptakeLV     Rate of N uptake in leaves                        Y   |kg N ha-1 d-1|
    RNuptakeST     Rate of N uptake in stems                         Y   |kg N ha-1 d-1|
    RNuptakeRT     Rate of N uptake in roots                         Y   |kg N ha-1 d-1|
    RNuptakeSO     Rate of N uptake in storage organs                Y   |kg N ha-1 d-1|

    RPuptakeLV     Rate of P uptake in leaves                        Y   |kg P ha-1 d-1|
    RPuptakeST     Rate of P uptake in stems                         Y   |kg P ha-1 d-1|
    RPuptakeRT     Rate of P uptake in roots                         Y   |kg P ha-1 d-1|
    RPuptakeSO     Rate of P uptake in storage organs                Y   |kg P ha-1 d-1|

    RKuptakeLV     Rate of K uptake in leaves                        Y   |kg K ha-1 d-1|
    RKuptakeST     Rate of K uptake in stems                         Y   |kg K ha-1 d-1|
    RKuptakeRT     Rate of K uptake in roots                         Y   |kg K ha-1 d-1|
    RKuptakeSO     Rate of K uptake in storage organs                Y   |kg K ha-1 d-1|

    RNuptake       Total rate of N uptake                            Y   |kg N ha-1 d-1|
    RPuptake       Total rate of P uptake                            Y   |kg P ha-1 d-1|
    RKuptake       Total rate of K uptake                            Y   |kg K ha-1 d-1|
    RNfixation     Rate of N fixation                                Y   |kg K ha-1 d-1|
    ===========  ================================================= ==== ================

    **Signals send or handled**

    None

    **External dependencies**

    ================  =================================== ====================  ===========
     Name              Description                         Provided by            Unit
    ================  =================================== ====================  ===========
    DVS               Crop development stage              DVS_Phenology              -
    TRA               Crop transpiration                  Evapotranspiration     |cm d-1|
    TRAMX             Potential crop transpiration        Evapotranspiration     |cm d-1|
    NAVAIL            Total available N from soil         NPK_Soil_Dynamics      |kg ha-1|
    PAVAIL            Total available P from soil         NPK_Soil_Dynamics      |kg ha-1|
    KAVAIL            Total available K from soil         NPK_Soil_Dynamics      |kg ha-1|
    Ntranslocatable   Translocatable amount of N from     NPK_Translocation      |kg ha-1|
                      stems, Leaves and roots
    Ptranslocatable   As for P                            NPK_Translocation      |kg ha-1|
    Ktranslocatable   As for K                            NPK_Translocation      |kg ha-1|
    ================  =================================== ====================  ===========

    """

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        PMAXLV_TB = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        KMAXLV_TB = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        
        NMAXRT_FR = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        PMAXRT_FR = Float(-99.)  # maximum P concentration in roots as fraction of maximum P concentration in leaves
        KMAXRT_FR = Float(-99.)  # maximum K concentration in roots as fraction of maximum K concentration in leaves

        NMAXST_FR = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves
        PMAXST_FR = Float(-99.)  # maximum P concentration in stems as fraction of maximum P concentration in leaves
        KMAXST_FR = Float(-99.)  # maximum K concentration in stems as fraction of maximum K concentration in leaves
        
        NMAXSO = Float(-99.)  # maximum P concentration in storage organs [kg N kg-1 dry biomass]
        PMAXSO = Float(-99.)  # maximum P concentration in storage organs [kg P kg-1 dry biomass]
        KMAXSO = Float(-99.)  # maximum K concentration in storage organs [kg K kg-1 dry biomass]
        
        TCNT = Float(-99.)  # time coefficient for N translocation to storage organs [days]
        TCPT = Float(-99.)  # time coefficient for P translocation to storage organs [days]
        TCKT = Float(-99.)  # time coefficient for K translocation to storage organs [days]

        NFIX_FR = Float(-99.)  # fraction of crop nitrogen uptake by biological fixation
        DVS_NPK_STOP = Float(-99.)  # development stage above which no crop N-P-K uptake does occur

    class RateVariables(RatesTemplate):
        RNuptakeLV = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RNuptakeST = Float(-99.)
        RNuptakeRT = Float(-99.)
        RNuptakeSO = Float(-99.)

        RPuptakeLV = Float(-99.)  # P uptake rate [kg ha-1 d -1]
        RPuptakeST = Float(-99.)
        RPuptakeRT = Float(-99.)
        RPuptakeSO = Float(-99.)

        RKuptakeLV = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RKuptakeST = Float(-99.)
        RKuptakeRT = Float(-99.)
        RKuptakeSO = Float(-99.)

        RNuptake = Float(-99.)  # Total N uptake rate [kg ha-1 d -1]
        RPuptake = Float(-99.)
        RKuptake = Float(-99.)
        RNfixation = Float(-99.)

        NdemandLV = Float(-99.)
        NdemandST = Float(-99.)
        NdemandRT = Float(-99.)
        NdemandSO = Float(-99.)

        PdemandLV = Float(-99.)
        PdemandST = Float(-99.)
        PdemandRT = Float(-99.)
        PdemandSO = Float(-99.)

        KdemandLV = Float(-99.)
        KdemandST = Float(-99.)
        KdemandRT = Float(-99.)
        KdemandSO = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: a ParameterProvider with parameter key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk

        self.rates = self.RateVariables(kiosk,
            publish=["RNuptakeLV", "RNuptakeST", "RNuptakeRT", "RNuptakeSO",
                     "RPuptakeLV", "RPuptakeST", "RPuptakeRT", "RPuptakeSO",
                     "RKuptakeLV", "RKuptakeST", "RKuptakeRT", "RKuptakeSO",
                     "RNuptake", "RPuptake", "RKuptake", "RNfixation"])

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        delt = 1.0
        mc = self._compute_NPK_max_concentrations()

        # Total NPK demand of leaves, stems, roots and storage organs
        # Demand consists of a demand carried over from previous timesteps plus a demand from new growth
        # Note that we are pre-integrating here, so a multiplication with time-step delt is required

        # N demand [kg ha-1]
        r.NdemandLV = max(mc.NMAXLV * k.WLV - k.NamountLV, 0.) + max(k.GRLV * mc.NMAXLV, 0) * delt
        r.NdemandST = max(mc.NMAXST * k.WST - k.NamountST, 0.) + max(k.GRST * mc.NMAXST, 0) * delt
        r.NdemandRT = max(mc.NMAXRT * k.WRT - k.NamountRT, 0.) + max(k.GRRT * mc.NMAXRT, 0) * delt
        r.NdemandSO = max(mc.NMAXSO * k.WSO - k.NamountSO, 0.)

        # P demand [kg ha-1]
        r.PdemandLV = max(mc.PMAXLV * k.WLV - k.PamountLV, 0.) + max(k.GRLV * mc.PMAXLV, 0) * delt
        r.PdemandST = max(mc.PMAXST * k.WST - k.PamountST, 0.) + max(k.GRST * mc.PMAXST, 0) * delt
        r.PdemandRT = max(mc.PMAXRT * k.WRT - k.PamountRT, 0.) + max(k.GRRT * mc.PMAXRT, 0) * delt
        r.PdemandSO = max(mc.PMAXSO * k.WSO - k.PamountSO, 0.)

        # K demand [kg ha-1]
        r.KdemandLV = max(mc.KMAXLV * k.WLV - k.KamountLV, 0.) + max(k.GRLV * mc.KMAXLV, 0) * delt
        r.KdemandST = max(mc.KMAXST * k.WST - k.KamountST, 0.) + max(k.GRST * mc.KMAXST, 0) * delt
        r.KdemandRT = max(mc.KMAXRT * k.WRT - k.KamountRT, 0.) + max(k.GRRT * mc.KMAXRT, 0) * delt
        r.KdemandSO = max(mc.KMAXSO * k.WSO - k.KamountSO, 0.)

        Ndemand = r.NdemandLV + r.NdemandST + r.NdemandRT
        Pdemand = r.PdemandLV + r.PdemandST + r.PdemandRT
        Kdemand = r.KdemandLV + r.KdemandST + r.KdemandRT

        # NPK uptake rate in storage organs (kg N ha-1 d-1) is the mimimum of supply and
        # demand divided by the time coefficient for N/P/K translocation
        r.RNuptakeSO = min(r.NdemandSO, k.Ntranslocatable)/p.TCNT
        r.RPuptakeSO = min(r.PdemandSO, k.Ptranslocatable)/p.TCPT
        r.RKuptakeSO = min(r.KdemandSO, k.Ktranslocatable)/p.TCKT

        # No nutrients are absorbed after development stage DVS_NPK_STOP or
        # when severe water shortage occurs i.e. RFTRA <= 0.01
        if k.DVS < p.DVS_NPK_STOP and k.RFTRA > 0.01:
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.

        # biological nitrogen fixation
        r.RNfixation = (max(0., p.NFIX_FR * Ndemand) * NutrientLIMIT)

        # NPK uptake rate from soil
        r.RNuptake = (max(0., min(Ndemand - r.RNfixation, k.NAVAIL)) * NutrientLIMIT)
        r.RPuptake = (max(0., min(Pdemand, k.PAVAIL)) * NutrientLIMIT)
        r.RKuptake = (max(0., min(Kdemand, k.KAVAIL)) * NutrientLIMIT)

        # NPK uptake rate
        # if no demand then uptake rate = 0.
        if Ndemand == 0.:
            r.RNuptakeLV = r.RNuptakeST = r.RNuptakeRT = 0.
        else:
            r.RNuptakeLV = (r.NdemandLV / Ndemand) * (r.RNuptake + r.RNfixation)
            r.RNuptakeST = (r.NdemandST / Ndemand) * (r.RNuptake + r.RNfixation)
            r.RNuptakeRT = (r.NdemandRT / Ndemand) * (r.RNuptake + r.RNfixation)

        if Pdemand == 0.:
            r.RPuptakeLV = r.RPuptakeST = r.RPuptakeRT = 0.
        else:
            r.RPuptakeLV = (r.PdemandLV / Pdemand) * r.RPuptake
            r.RPuptakeST = (r.PdemandST / Pdemand) * r.RPuptake
            r.RPuptakeRT = (r.PdemandRT / Pdemand) * r.RPuptake

        if Kdemand == 0.:
            r.RKuptakeLV = r.RKuptakeST = r.RKuptakeRT = 0.
        else:
            r.RKuptakeLV = (r.KdemandLV / Kdemand) * r.RKuptake
            r.RKuptakeST = (r.KdemandST / Kdemand) * r.RKuptake
            r.RKuptakeRT = (r.KdemandRT / Kdemand) * r.RKuptake

    @prepare_states
    def integrate(self, day, delt=1.0):
        pass

    def _compute_NPK_max_concentrations(self):
        """Computes the maximum N/P/K concentrations in leaves, stems, roots and storage organs.
        
        Note that max concentrations are first derived from the dilution curve for leaves. 
        Maximum concentrations for stems and roots are computed as a fraction of the 
        concentration for leaves. Maximum concentration for storage organs is directly taken from
        the parameters N/P/KMAXSO.
        """

        p = self.params
        k = self.kiosk
        NMAXLV = p.NMAXLV_TB(k.DVS)
        PMAXLV = p.PMAXLV_TB(k.DVS)
        KMAXLV = p.KMAXLV_TB(k.DVS)
        max_NPK_conc = MaxNutrientConcentrations(
            # Maximum NPK concentrations in leaves [kg N kg-1 DM]
            NMAXLV=NMAXLV,
            PMAXLV=PMAXLV,
            KMAXLV=KMAXLV,
            # Maximum NPK concentrations in stems and roots [kg N kg-1 DM]
            NMAXST=(p.NMAXST_FR * NMAXLV),
            NMAXRT=p.NMAXRT_FR * NMAXLV,
            NMAXSO=p.NMAXSO,

            PMAXST=p.PMAXST_FR * PMAXLV,
            PMAXRT=p.PMAXRT_FR * PMAXLV,
            PMAXSO=p.PMAXSO,

            KMAXST=p.KMAXST_FR * KMAXLV,
            KMAXRT=p.KMAXRT_FR * KMAXLV,
            KMAXSO=p.KMAXSO
        )

        return max_NPK_conc
