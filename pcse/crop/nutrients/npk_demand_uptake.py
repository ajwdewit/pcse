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

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    NDEMLV     N Demand in living leaves                         N   |kg N ha-1|
    NDEMST     N Demand in living stems                          N   |kg N ha-1|
    NDEMRT     N Demand in living roots                          N   |kg N ha-1|
    NDEMSO     N Demand in storage organs                        N   |kg N ha-1|

    PDEMLV     P Demand in living leaves                         N   |kg P ha-1|
    PDEMST     P Demand in living stems                          N   |kg P ha-1|
    PDEMRT     P Demand in living roots                          N   |kg P ha-1|
    PDEMSO     P Demand in storage organs                        N   |kg P ha-1|

    KDEMLV     K Demand in living leaves                         N   |kg K ha-1|
    KDEMST     K Demand in living stems                          N   |kg K ha-1|
    KDEMRT     K Demand in living roots                          N   |kg K ha-1|
    KDEMSO     K Demand in storage organs                        N   |kg K ha-1|
    =======  ================================================= ==== ============


    **Rate variables**

    =======  ================================================= ==== ================
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ================
    RNULV     Rate of N uptake in leaves                         Y   |kg N ha-1 d-1|
    RNUST     Rate of N uptake in stems                          Y   |kg N ha-1 d-1|
    RNURT     Rate of N uptake in roots                          Y   |kg N ha-1 d-1|
    RNUSO     Rate of N uptake in storage organs                 Y   |kg N ha-1 d-1|

    RPULV     Rate of P uptake in leaves                         Y   |kg P ha-1 d-1|
    RPUST     Rate of P uptake in stems                          Y   |kg P ha-1 d-1|
    RPURT     Rate of P uptake in roots                          Y   |kg P ha-1 d-1|
    RPUSO     Rate of P uptake in storage organs                 Y   |kg P ha-1 d-1|

    RKULV     Rate of K uptake in leaves                         Y   |kg K ha-1 d-1|
    RKUST     Rate of K uptake in stems                          Y   |kg K ha-1 d-1|
    RKURT     Rate of K uptake in roots                          Y   |kg K ha-1 d-1|
    RKUSO     Rate of K uptake in storage organs                 Y   |kg K ha-1 d-1|

    RNUPTAKE  Total rate of N uptake                             Y   |kg N ha-1 d-1|
    RPUPTAKE  Total rate of P uptake                             Y   |kg P ha-1 d-1|
    RKUPTAKE  Total rate of K uptake                             Y   |kg K ha-1 d-1|
    RNFIX     Rate of N fixation                                 Y   |kg K ha-1 d-1|
    =======  ================================================= ==== ================

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
    NTRANSLOCATABLE   Translocatable amount of N from     NPK_Translocation      |kg ha-1|
                      stems, Leaves and roots
    PTRANSLOCATABLE   As for P                            NPK_Translocation      |kg ha-1|
    KTRANSLOCATABLE   As for K                            NPK_Translocation      |kg ha-1|
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
        NDEMLV = Float(-99.)
        NDEMST = Float(-99.)
        NDEMRT = Float(-99.)
        NDEMSO = Float(-99.)

        PDEMLV = Float(-99.)
        PDEMST = Float(-99.)
        PDEMRT = Float(-99.)
        PDEMSO = Float(-99.)
        
        KDEMLV = Float(-99.)
        KDEMST = Float(-99.)
        KDEMRT = Float(-99.)
        KDEMSO = Float(-99.)

    class RateVariables(RatesTemplate):
        RNULV = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RNUST = Float(-99.)
        RNURT = Float(-99.)
        RNUSO = Float(-99.)

        RPULV = Float(-99.)  # P uptake rate [kg ha-1 d -1]
        RPUST = Float(-99.)
        RPURT = Float(-99.)
        RPUSO = Float(-99.)

        RKULV = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RKUST = Float(-99.)
        RKURT = Float(-99.)
        RKUSO = Float(-99.)

        RNUPTAKE = Float(-99.)  # Total N uptake rate [kg ha-1 d -1]
        RPUPTAKE = Float(-99.)
        RKUPTAKE = Float(-99.)
        RNFIX = Float(-99.)

    
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: a ParameterProvider with parameter key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk

        self.rates = self.RateVariables(kiosk,
            publish=["RNULV", "RNUST", "RNURT", "RNUSO",
                     "RPULV", "RPUST", "RPURT", "RPUSO",
                     "RKULV", "RKUST", "RKURT", "RKUSO",
                     "RNUPTAKE", "RPUPTAKE", "RKUPTAKE", "RNFIX"])

        self.states = self.StateVariables(kiosk,
            NDEMLV=0., NDEMST=0., NDEMRT=0., NDEMSO=0.,
            PDEMLV=0., PDEMST=0., PDEMRT=0., PDEMSO=0.,
            KDEMLV=0., KDEMST=0., KDEMRT=0., KDEMSO=0.)

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

#       total NPK demand of leaves, stems and roots
        NDEMTO = s.NDEMLV + s.NDEMST + s.NDEMRT
        PDEMTO = s.PDEMLV + s.PDEMST + s.PDEMRT
        KDEMTO = s.KDEMLV + s.KDEMST + s.KDEMRT

#       NPK uptake rate in storage organs (kg N ha-1 d-1)
#       is the mimimum of supply and demand divided by the
#       time coefficient for N/P/K translocation
        r.RNUSO = min(s.NDEMSO, k.NTRANSLOCATABLE)/p.TCNT
        r.RPUSO = min(s.PDEMSO, k.PTRANSLOCATABLE)/p.TCPT
        r.RKUSO = min(s.KDEMSO, k.KTRANSLOCATABLE)/p.TCKT

#       No nutrients are absorbed after development stage DVS_NPK_STOP or
#       when severe water shortage occurs i.e. RFTRA <= 0.01
        if k.DVS < p.DVS_NPK_STOP and k.RFTRA > 0.01:
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.

        # biological nitrogen fixation
        r.RNFIX = (max(0., p.NFIX_FR * NDEMTO) * NutrientLIMIT)

        # NPK uptake rate from soil
        r.RNUPTAKE = (max(0., min(NDEMTO - r.RNFIX, k.NAVAIL)) * NutrientLIMIT)
        r.RPUPTAKE = (max(0., min(PDEMTO, k.PAVAIL)) * NutrientLIMIT)
        r.RKUPTAKE = (max(0., min(KDEMTO, k.KAVAIL)) * NutrientLIMIT)

        # NPK uptake rate
        # if no demand then uptake rate = 0.
        if NDEMTO == 0.:
            r.RNULV = r.RNUST = r.RNURT = 0.
        else:
            r.RNULV = (s.NDEMLV / NDEMTO) * (r.RNUPTAKE + r.RNFIX)
            r.RNUST = (s.NDEMST / NDEMTO) * (r.RNUPTAKE + r.RNFIX)
            r.RNURT = (s.NDEMRT / NDEMTO) * (r.RNUPTAKE + r.RNFIX)

        if PDEMTO == 0.:
            r.RPULV = r.RPUST = r.RPURT = 0.
        else:
            r.RPULV = (s.PDEMLV / PDEMTO) * r.RPUPTAKE
            r.RPUST = (s.PDEMST / PDEMTO) * r.RPUPTAKE
            r.RPURT = (s.PDEMRT / PDEMTO) * r.RPUPTAKE

        if KDEMTO == 0.:
            r.RKULV = r.RKUST = r.RKURT = 0.
        else:
            r.RKULV = (s.KDEMLV / KDEMTO) * r.RKUPTAKE
            r.RKUST = (s.KDEMST / KDEMTO) * r.RKUPTAKE
            r.RKURT = (s.KDEMRT / KDEMTO) * r.RKUPTAKE

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
        s.NDEMLV = max(NMAXLV * k.WLV - k.ANLV, 0.)  # maybe should be divided by one day, see equation 5 Shibu etal 2010
        s.NDEMST = max(NMAXST * k.WST - k.ANST, 0.)
        s.NDEMRT = max(NMAXRT * k.WRT - k.ANRT, 0.)
        s.NDEMSO = max(NMAXSO * k.WSO - k.ANSO, 0.)

#       P demand [kg ha-1]
        s.PDEMLV = max(PMAXLV * k.WLV - k.APLV, 0.)
        s.PDEMST = max(PMAXST * k.WST - k.APST, 0.)
        s.PDEMRT = max(PMAXRT * k.WRT - k.APRT, 0.)
        s.PDEMSO = max(PMAXSO * k.WSO - k.APSO, 0.)

#       K demand [kg ha-1]
        s.KDEMLV = max(KMAXLV * k.WLV - k.AKLV, 0.)
        s.KDEMST = max(KMAXST * k.WST - k.AKST, 0.)
        s.KDEMRT = max(KMAXRT * k.WRT - k.AKRT, 0.)
        s.KDEMSO = max(KMAXSO * k.WSO - k.AKSO, 0.)
