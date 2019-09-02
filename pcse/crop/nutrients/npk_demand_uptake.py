# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf

from ...base import StatesTemplate, ParamTemplate, SimulationObject, RatesTemplate
from ...decorators import prepare_rates, prepare_states
from ...traitlets import HasTraits, Float, Int, Instance
from ...util import AfgenTrait

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

    ============  ============================================= =======  ======================
     Name          Description                                   Type     Unit
    ============  ============================================= =======  ======================
    NMAXLV_TB      Maximum N concentration in leaves as          TCr     kg N kg-1 dry biomass
                   function of DVS
    PMAXLV_TB      As for P                                      TCr     kg P kg-1 dry biomass
    KMAXLV_TB      As for K                                      TCr     kg K kg-1 dry biomass

    NMAXRT_FR      Maximum N concentration in roots as fraction  SCr     -
                   of maximum N concentration in leaves
    PMAXRT_FR      As for P                                      SCr     -
    KMAXRT_FR      As for K                                      SCr     -

    NMAXST_FR      Maximum N concentration in stems as fraction  SCr     -
                   of maximum N concentration in leaves
    PMAXST_FR      As for P                                      SCr     -
    KMAXST_FR      As for K                                      SCr     -

    NMAXSO         Maximum N concentration in storage organs     SCr     kg N kg-1 dry biomass
    PMAXSO         As for P                                      SCr     kg P kg-1 dry biomass
    KMAXSO         As for K                                      SCr     kg K kg-1 dry biomass

    NCRIT_FR       Critical N concentration as fraction of       SCr     -
                   maximum N concentration for vegetative
                   plant organs as a whole (leaves + stems)
    PCRIT_FR       As for P                                      SCr     -
    KCRIT_FR       As for K                                      SCr     -

    TCNT           Time coefficient for N translation to         SCr     days
                   storage organs
    TCPT           As for P                                      SCr     days
    TCKT           As for K                                      SCr     days

    NFIX_FR        fraction of crop nitrogen uptake by           SCr     kg N kg-1 dry biomass
                   biological fixation
    DVS_NPK_STOP    Development stage after which no nutrients    SCr     -
                   are taken up from the soil by the crop.
    ============  ============================================= ======= =======================

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
    RNuptakeLV     Rate of N uptake in leaves                         Y   |kg N ha-1 d-1|
    RNuptakeST     Rate of N uptake in stems                          Y   |kg N ha-1 d-1|
    RNuptakeRT     Rate of N uptake in roots                          Y   |kg N ha-1 d-1|
    RNuptakeSO     Rate of N uptake in storage organs                 Y   |kg N ha-1 d-1|

    RPuptakeLV     Rate of P uptake in leaves                         Y   |kg P ha-1 d-1|
    RPuptakeST     Rate of P uptake in stems                          Y   |kg P ha-1 d-1|
    RPuptakeRT     Rate of P uptake in roots                          Y   |kg P ha-1 d-1|
    RPuptakeSO     Rate of P uptake in storage organs                 Y   |kg P ha-1 d-1|

    RKuptakeLV     Rate of K uptake in leaves                         Y   |kg K ha-1 d-1|
    RKuptakeST     Rate of K uptake in stems                          Y   |kg K ha-1 d-1|
    RKuptakeRT     Rate of K uptake in roots                          Y   |kg K ha-1 d-1|
    RKuptakeSO     Rate of K uptake in storage organs                 Y   |kg K ha-1 d-1|

    RNuptake       Total rate of N uptake                             Y   |kg N ha-1 d-1|
    RPuptake       Total rate of P uptake                             Y   |kg P ha-1 d-1|
    RKuptake       Total rate of K uptake                             Y   |kg K ha-1 d-1|
    RNfixation     Rate of N fixation                                 Y   |kg K ha-1 d-1|
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

    class StateVariables(StatesTemplate):
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

        self.states = self.StateVariables(kiosk,
            NdemandLV=0., NdemandST=0., NdemandRT=0., NdemandSO=0.,
            PdemandLV=0., PdemandST=0., PdemandRT=0., PdemandSO=0.,
            KdemandLV=0., KdemandST=0., KdemandRT=0., KdemandSO=0.)

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        # total NPK demand of leaves, stems and roots
        Ndemand = s.NdemandLV + s.NdemandST + s.NdemandRT
        Pdemand = s.PdemandLV + s.PdemandST + s.PdemandRT
        Kdemand = s.KdemandLV + s.KdemandST + s.KdemandRT

        # NPK uptake rate in storage organs (kg N ha-1 d-1) is the mimimum of supply and
        # demand divided by the time coefficient for N/P/K translocation
        r.RNuptakeSO = min(s.NdemandSO, k.Ntranslocatable)/p.TCNT
        r.RPuptakeSO = min(s.PdemandSO, k.Ptranslocatable)/p.TCPT
        r.RKuptakeSO = min(s.KdemandSO, k.Ktranslocatable)/p.TCKT

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
            r.RNuptakeLV = (s.NdemandLV / Ndemand) * (r.RNuptake + r.RNfixation)
            r.RNuptakeST = (s.NdemandST / Ndemand) * (r.RNuptake + r.RNfixation)
            r.RNuptakeRT = (s.NdemandRT / Ndemand) * (r.RNuptake + r.RNfixation)

        if Pdemand == 0.:
            r.RPuptakeLV = r.RPuptakeST = r.RPuptakeRT = 0.
        else:
            r.RPuptakeLV = (s.PdemandLV / Pdemand) * r.RPuptake
            r.RPuptakeST = (s.PdemandST / Pdemand) * r.RPuptake
            r.RPuptakeRT = (s.PdemandRT / Pdemand) * r.RPuptake

        if Kdemand == 0.:
            r.RKuptakeLV = r.RKuptakeST = r.RKuptakeRT = 0.
        else:
            r.RKuptakeLV = (s.KdemandLV / Kdemand) * r.RKuptake
            r.RKuptakeST = (s.KdemandST / Kdemand) * r.RKuptake
            r.RKuptakeRT = (s.KdemandRT / Kdemand) * r.RKuptake

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        p = self.params
        k = self.kiosk

#       Maximum NPK concentrations in leaves [kg N kg-1 DM]
        NMAXLV = p.NMAXLV_TB(k.DVS)
        PMAXLV = p.PMAXLV_TB(k.DVS)
        KMAXLV = p.KMAXLV_TB(k.DVS)
        
#       Maximum NPK concentrations in stems and roots [kg N kg-1 DM]
        NMAXST = p.NMAXST_FR * NMAXLV
        NMAXRT = p.NMAXRT_FR * NMAXLV
        NMAXSO = p.NMAXSO
      
        PMAXST = p.PMAXST_FR * PMAXLV
        PMAXRT = p.PMAXRT_FR * PMAXLV
        PMAXSO = p.PMAXSO
      
        KMAXST = p.KMAXST_FR * KMAXLV
        KMAXRT = p.KMAXRT_FR * KMAXLV
        KMAXSO = p.KMAXSO

#       N demand [kg ha-1]
        s.NdemandLV = max(NMAXLV * k.WLV - k.NamountLV, 0.)  # maybe should be divided by one day, see equation 5 Shibu etal 2010
        s.NdemandST = max(NMAXST * k.WST - k.NamountST, 0.)
        s.NdemandRT = max(NMAXRT * k.WRT - k.NamountRT, 0.)
        s.NdemandSO = max(NMAXSO * k.WSO - k.NamountSO, 0.)

#       P demand [kg ha-1]
        s.PdemandLV = max(PMAXLV * k.WLV - k.PamountLV, 0.)
        s.PdemandST = max(PMAXST * k.WST - k.PamountST, 0.)
        s.PdemandRT = max(PMAXRT * k.WRT - k.PamountRT, 0.)
        s.PdemandSO = max(PMAXSO * k.WSO - k.PamountSO, 0.)

#       K demand [kg ha-1]
        s.KdemandLV = max(KMAXLV * k.WLV - k.KamountLV, 0.)
        s.KdemandST = max(KMAXST * k.WST - k.KamountST, 0.)
        s.KdemandRT = max(KMAXRT * k.WRT - k.KamountRT, 0.)
        s.KdemandSO = max(KMAXSO * k.WSO - k.KamountSO, 0.)
