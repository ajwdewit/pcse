#!/usr/bin/env python

from .. import exceptions as exc
from ..traitlets import Float, Int, Instance
from ..decorators import prepare_rates, prepare_states
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from ..util import AfgenTrait
from .nutrients import NPK_Translocation
from .nutrients import NPK_Demand_Uptake

class NPK_Crop_Dynamics(SimulationObject):
    """Implementation of overall NPK crop dynamics.

    NPK_Crop_Dynamics implements the overall logic of N/P/K book-keeping within the
    crop.

    **Simulation parameters**
    
    =============  ================================================= =======================
     Name           Description                                        Unit
    =============  ================================================= =======================
    DVS_NPK_STOP   DVS above which no crop N-P-K uptake occurs           -

    NMAXLV_TB      Maximum N concentration in leaves as               kg N kg-1 dry biomass
                   function of dvs
    PMAXLV_TB      As for P                                           kg P kg-1 dry biomass
    KMAXLV_TB      As for K                                           kg K kg-1 dry biomass

    NMAXRT_FR      Maximum N concentration in roots as fraction       -
                   of maximum N concentration in leaves
    PMAXRT_FR      As for P                                           -
    KMAXRT_FR      As for K                                           -

    NMAXST_FR      Maximum N concentration in stems as fraction       -
                   of maximum N concentration in leaves
    KMAXST_FR      As for K                                           -
    PMAXST_FR      As for P                                           -

    NRESIDLV       Residual N fraction in leaves                      kg N kg-1 dry biomass
    PRESIDLV       Residual P fraction in leaves                      kg P kg-1 dry biomass
    KRESIDLV       Residual K fraction in leaves                      kg K kg-1 dry biomass

    NRESIDRT       Residual N fraction in roots                       kg N kg-1 dry biomass
    PRESIDRT       Residual P fraction in roots                       kg P kg-1 dry biomass
    KRESIDRT       Residual K fraction in roots                       kg K kg-1 dry biomass

    NRESIDST       Residual N fraction in stems                       kg N kg-1 dry biomass
    PRESIDST       Residual P fraction in stems                       kg P kg-1 dry biomass
    KRESIDST       Residual K fraction in stems                       kg K kg-1 dry biomass
    =============  ================================================= =======================

    **State variables**

    ==========  ================================================== ============
     Name        Description                                          Unit
    ==========  ================================================== ============
    NamountLV     Actual N amount in living leaves                  |kg N ha-1|
    PamountLV     Actual P amount in living leaves                  |kg P ha-1|
    KamountLV     Actual K amount in living leaves                  |kg K ha-1|
        
    NamountST     Actual N amount in living stems                   |kg N ha-1|
    PamountST     Actual P amount in living stems                   |kg P ha-1|
    KamountST     Actual K amount in living stems                   |kg K ha-1|

    NamountSO     Actual N amount in living storage organs          |kg N ha-1|
    PamountSO     Actual P amount in living storage organs          |kg P ha-1|
    KamountSO     Actual K amount in living storage organs          |kg K ha-1|
    
    NamountRT     Actual N amount in living roots                   |kg N ha-1|
    PamountRT     Actual P amount in living roots                   |kg P ha-1|
    KamountRT     Actual K amount in living roots                   |kg K ha-1|
    
    Nuptake_T    total absorbed N amount                            |kg N ha-1|
    Puptake_T    total absorbed P amount                            |kg P ha-1|
    Kuptake_T    total absorbed K amount                            |kg K ha-1|
    Nfix_T       total biological fixated N amount                  |kg N ha-1|
    ==========  ================================================== ============

    **Rate variables**

    ===========  =================================================  ================
     Name         Description                                           Unit
    ===========  =================================================  ================
    RNamountLV     Weight increase (N) in leaves                    |kg N ha-1 d-1|
    RPamountLV     Weight increase (P) in leaves                    |kg P ha-1 d-1|
    RKamountLV     Weight increase (K) in leaves                    |kg K ha-1 d-1|
    
    RNamountST     Weight increase (N) in stems                     |kg N ha-1 d-1|
    RPamountST     Weight increase (P) in stems                     |kg P ha-1 d-1|
    RKamountST     Weight increase (K) in stems                     |kg K ha-1 d-1|
        
    RNamountRT     Weight increase (N) in roots                     |kg N ha-1 d-1|
    RPamountRT     Weight increase (P) in roots                     |kg P ha-1 d-1|
    RKamountRT     Weight increase (K) in roots                     |kg K ha-1 d-1|
    
    RNamountSO     Weight increase (N) in storage organs            |kg N ha-1 d-1|
    RPamountSO     Weight increase (P) in storage organs            |kg P ha-1 d-1|
    RKamountSO     Weight increase (K) in storage organs            |kg K ha-1 d-1|

    RNdeathLV      Rate of N loss in leaves                         |kg N ha-1 d-1|
    RPdeathLV      as for P                                         |kg P ha-1 d-1|
    RKdeathLV      as for K                                         |kg K ha-1 d-1|

    RNdeathST      Rate of N loss in roots                          |kg N ha-1 d-1|
    RPdeathST      as for P                                         |kg P ha-1 d-1|
    RKdeathST      as for K                                         |kg K ha-1 d-1|

    RNdeathRT      Rate of N loss in stems                          |kg N ha-1 d-1|
    RPdeathRT      as for P                                         |kg P ha-1 d-1|
    RKdeathRT      as for K                                         |kg K ha-1 d-1|

    RNloss         N loss due to senescence                         |kg N ha-1 d-1|
    RPloss         P loss due to senescence                         |kg P ha-1 d-1|
    RKloss         K loss due to senescence                         |kg K ha-1 d-1|
    ===========  =================================================  ================
    
    **Signals send or handled**
    
    None
    
    **External dependencies**
    
    =======  =================================== ====================  ==============
     Name     Description                         Provided by            Unit
    =======  =================================== ====================  ==============
    DVS      Crop development stage              DVS_Phenology           -
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    WRT      Dry weight of living roots          WOFOST_Root_Dynamics  |kg ha-1|
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    DRLV     Death rate of leaves                WOFOST_Leaf_Dynamics  |kg ha-1 d-1|
    DRRT     Death rate of roots                 WOFOST_Root_Dynamics  |kg ha-1 d-1|
    DRST     Death rate of stems                 WOFOST_Stem_Dynamics  |kg ha-1 d-1|
    =======  =================================== ====================  ==============
    """

    translocation = Instance(SimulationObject)
    demand_uptake = Instance(SimulationObject)

    NamountLVI = Float(-99.)  # initial soil N amount in leaves
    NamountSTI = Float(-99.)  # initial soil N amount in stems
    NamountRTI = Float(-99.)  # initial soil N amount in roots
    NamountSOI = Float(-99.)  # initial soil N amount in storage organs
    
    PamountLVI = Float(-99.)  # initial soil P amount in leaves
    PamountSTI = Float(-99.)  # initial soil P amount in stems
    PamountRTI = Float(-99.)  # initial soil P amount in roots
    PamountSOI = Float(-99.)  # initial soil P amount in storage organs

    KamountLVI = Float(-99.)  # initial soil K amount in leaves
    KamountSTI = Float(-99.)  # initial soil K amount in stems
    KamountRTI = Float(-99.)  # initial soil K amount in roots
    KamountSOI = Float(-99.)  # initial soil K amount in storage organs

    class Parameters(ParamTemplate):
        DVS_NPK_STOP = Float(-99.)
        NMAXLV_TB = AfgenTrait()
        PMAXLV_TB = AfgenTrait()
        KMAXLV_TB = AfgenTrait()
        NMAXST_FR = Float(-99.)
        NMAXRT_FR = Float(-99.)
        PMAXST_FR = Float(-99.)
        PMAXRT_FR = Float(-99.)
        KMAXST_FR = Float(-99.)
        KMAXRT_FR = Float(-99.)
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NRESIDRT = Float(-99.)  # residual N fraction in roots [kg N kg-1 dry biomass]
        PRESIDLV = Float(-99.)  # residual P fraction in leaves [kg P kg-1 dry biomass]
        PRESIDST = Float(-99.)  # residual P fraction in stems [kg P kg-1 dry biomass]
        PRESIDRT = Float(-99.)  # residual P fraction in roots [kg P kg-1 dry biomass]
        KRESIDLV = Float(-99.)  # residual K fraction in leaves [kg K kg-1 dry biomass]
        KRESIDST = Float(-99.)  # residual K fraction in stems [kg K kg-1 dry biomass]
        KRESIDRT = Float(-99.)  # residual K fraction in roots [kg K kg-1 dry biomass]

    class StateVariables(StatesTemplate):
        NamountLV = Float(-99.) # N amount in leaves [kg N ha-1]
        PamountLV = Float(-99.) # P amount in leaves [kg P ]
        KamountLV = Float(-99.) # K amount in leaves [kg K ]
        
        NamountST = Float(-99.) # N amount in stems [kg N ]
        PamountST = Float(-99.) # P amount in stems [kg P ]
        KamountST = Float(-99.) # K amount in stems [kg K ]
      
        NamountSO = Float(-99.) # N amount in storage organs [kg N ]
        PamountSO = Float(-99.) # P amount in storage organs [kg P ]
        KamountSO = Float(-99.) # K amount in storage organs [kg K ]
        
        NamountRT = Float(-99.) # N amount in roots [kg N ]
        PamountRT = Float(-99.) # P amount in roots [kg P ]
        KamountRT = Float(-99.) # K amount in roots [kg K ]
        
        Nuptake_T = Float(-99.) # total absorbed N amount [kg N ]
        Puptake_T = Float(-99.) # total absorbed P amount [kg P ]
        Kuptake_T = Float(-99.) # total absorbed K amount [kg K ]
        Nfix_T = Float(-99.) # total biological fixated N amount [kg N ]
        
        Nlosses_T = Float(-99.)
        Plosses_T = Float(-99.)
        Klosses_T = Float(-99.)

    class RateVariables(RatesTemplate):
        RNamountLV = Float(-99.)  # Net rates of NPK in different plant organs 
        RPamountLV = Float(-99.)
        RKamountLV = Float(-99.)
        
        RNamountST = Float(-99.)
        RPamountST = Float(-99.)
        RKamountST = Float(-99.)
               
        RNamountRT = Float(-99.)
        RPamountRT = Float(-99.)
        RKamountRT = Float(-99.)
        
        RNamountSO = Float(-99.)
        RPamountSO = Float(-99.)
        RKamountSO = Float(-99.)
               
        RNdeathLV = Float(-99.)  # N loss rate leaves [kg ha-1 d-1]
        RNdeathST = Float(-99.)  # N loss rate stems  [kg ha-1 d-1]
        RNdeathRT = Float(-99.)  # N loss rate roots  [kg ha-1 d-1]
        
        RPdeathLV = Float(-99.)  # P loss rate leaves [kg ha-1 d-1]
        RPdeathST = Float(-99.)  # P loss rate stems  [kg ha-1 d-1]
        RPdeathRT = Float(-99.)  # P loss rate roots  [kg ha-1 d-1]
        
        RKdeathLV = Float(-99.)  # K loss rate leaves [kg ha-1 d-1]
        RKdeathST = Float(-99.)  # K loss rate stems  [kg ha-1 d-1]
        RKdeathRT = Float(-99.)  # K loss rate roots  [kg ha-1 d-1]

        RNloss = Float(-99.)
        RPloss = Float(-99.)
        RKloss = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameters as key/value pairs
        """  
        
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # Initialize components of the npk_crop_dynamics
        self.translocation = NPK_Translocation(day, kiosk, parvalues)
        self.demand_uptake = NPK_Demand_Uptake(day, kiosk, parvalues)

        # INITIAL STATES
        params = self.params
        k = kiosk

        # Initial amounts
        self.NamountLVI = NamountLV = k.WLV * params.NMAXLV_TB(k.DVS)
        self.NamountSTI = NamountST = k.WST * params.NMAXLV_TB(k.DVS) * params.NMAXST_FR
        self.NamountRTI = NamountRT = k.WRT * params.NMAXLV_TB(k.DVS) * params.NMAXRT_FR
        self.NamountSOI = NamountSO = 0.
        
        self.PamountLVI = PamountLV = k.WLV * params.PMAXLV_TB(k.DVS)
        self.PamountSTI = PamountST = k.WST * params.PMAXLV_TB(k.DVS) * params.PMAXST_FR
        self.PamountRTI = PamountRT = k.WRT * params.PMAXLV_TB(k.DVS) * params.PMAXRT_FR
        self.PamountSOI = PamountSO = 0.

        self.KamountLVI = KamountLV = k.WLV * params.KMAXLV_TB(k.DVS)
        self.KamountSTI = KamountST = k.WST * params.KMAXLV_TB(k.DVS) * params.KMAXST_FR
        self.KamountRTI = KamountRT = k.WRT * params.KMAXLV_TB(k.DVS) * params.KMAXRT_FR
        self.KamountSOI = KamountSO = 0.

        self.states = self.StateVariables(kiosk,
                        publish=["NamountLV", "NamountST", "NamountRT", "NamountSO", "PamountLV", "PamountST",
                                 "PamountRT", "PamountSO", "KamountLV", "KamountST", "KamountRT", "KamountSO"],
                        NamountLV=NamountLV, NamountST=NamountST, NamountRT=NamountRT, NamountSO=NamountSO,
                        PamountLV=PamountLV, PamountST=PamountST, PamountRT=PamountRT, PamountSO=PamountSO,
                        KamountLV=KamountLV, KamountST=KamountST, KamountRT=KamountRT, KamountSO=KamountSO,
                        Nuptake_T=0, Puptake_T=0., Kuptake_T=0., Nfix_T=0.,
                        Nlosses_T=0, Plosses_T=0., Klosses_T=0.)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        params = self.params
        k = self.kiosk
        
        self.demand_uptake.calc_rates(day, drv)
        self.translocation.calc_rates(day, drv)

        # Compute loss of NPK due to death of plant material
        rates.RNdeathLV = params.NRESIDLV * k.DRLV
        rates.RNdeathST = params.NRESIDST * k.DRST
        rates.RNdeathRT = params.NRESIDRT * k.DRRT

        rates.RPdeathLV = params.PRESIDLV * k.DRLV
        rates.RPdeathST = params.PRESIDST * k.DRST
        rates.RPdeathRT = params.PRESIDRT * k.DRRT

        rates.RKdeathLV = params.KRESIDLV * k.DRLV
        rates.RKdeathST = params.KRESIDST * k.DRST
        rates.RKdeathRT = params.KRESIDRT * k.DRRT

        # N rates in leaves, stems, root and storage organs computed as
        # uptake - translocation - death.
        # except for storage organs which only take up as a result of translocation.
        rates.RNamountLV = k.RNuptakeLV - k.RNtranslocationLV - rates.RNdeathLV
        rates.RNamountST = k.RNuptakeST - k.RNtranslocationST - rates.RNdeathST
        rates.RNamountRT = k.RNuptakeRT - k.RNtranslocationRT - rates.RNdeathRT
        rates.RNamountSO = k.RNuptakeSO
        
        # P rates in leaves, stems, root and storage organs
        rates.RPamountLV = k.RPuptakeLV - k.RPtranslocationLV - rates.RPdeathLV
        rates.RPamountST = k.RPuptakeST - k.RPtranslocationST - rates.RPdeathST
        rates.RPamountRT = k.RPuptakeRT - k.RPtranslocationRT - rates.RPdeathRT
        rates.RPamountSO = k.RPuptakeSO

        # K rates in leaves, stems, root and storage organs
        rates.RKamountLV = k.RKuptakeLV - k.RKtranslocationLV - rates.RKdeathLV
        rates.RKamountST = k.RKuptakeST - k.RKtranslocationST - rates.RKdeathST
        rates.RKamountRT = k.RKuptakeRT - k.RKtranslocationRT - rates.RKdeathRT
        rates.RKamountSO = k.RKuptakeSO
        
        rates.RNloss = rates.RNdeathLV + rates.RNdeathST + rates.RNdeathRT
        rates.RPloss = rates.RPdeathLV + rates.RPdeathST + rates.RPdeathRT
        rates.RKloss = rates.RKdeathLV + rates.RKdeathST + rates.RKdeathRT

        self._check_N_balance(day)
        self._check_P_balance(day)
        self._check_K_balance(day)
        
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states
        k = self.kiosk

        # N amount in leaves, stems, root and storage organs
        states.NamountLV += rates.RNamountLV
        states.NamountST += rates.RNamountST
        states.NamountRT += rates.RNamountRT
        states.NamountSO += rates.RNamountSO
        
        # P amount in leaves, stems, root and storage organs
        states.PamountLV += rates.RPamountLV
        states.PamountST += rates.RPamountST
        states.PamountRT += rates.RPamountRT
        states.PamountSO += rates.RPamountSO

        # K amount in leaves, stems, root and storage organs
        states.KamountLV += rates.RKamountLV
        states.KamountST += rates.RKamountST
        states.KamountRT += rates.RKamountRT
        states.KamountSO += rates.RKamountSO
        
        self.translocation.integrate(day, delt)
        self.demand_uptake.integrate(day, delt)

        # total NPK uptake from soil
        states.Nuptake_T += k.RNuptake
        states.Puptake_T += k.RPuptake
        states.Kuptake_T += k.RKuptake
        states.Nfix_T += k.RNfixation
        
        states.Nlosses_T += rates.RNloss
        states.Plosses_T += rates.RPloss
        states.Klosses_T += rates.RKloss

    def _check_N_balance(self, day):
        s = self.states
        checksum = abs(s.Nuptake_T + s.Nfix_T +
                       (self.NamountLVI + self.NamountSTI + self.NamountRTI + self.NamountSOI) -
                       (s.NamountLV + s.NamountST + s.NamountRT + s.NamountSO + s.Nlosses_T))

        if abs(checksum) >= 1.0:
            msg = "N flows not balanced on day %s\n" % day
            msg += "Checksum: %f, Nuptake_T: %f, Nfix_T: %f\n" % (checksum, s.Nuptake_T, s.Nfix_T)
            msg += "NamountLVI: %f, NamountSTI: %f, NamountRTI: %f, NamountSOI: %f\n"  % \
                   (self.NamountLVI, self.NamountSTI, self.NamountRTI, self.NamountSOI)
            msg += "NamountLV: %f, NamountST: %f, NamountRT: %f, NamountSO: %f\n" % \
                   (s.NamountLV, s.NamountST, s.NamountRT, s.NamountSO)
            msg += "NLOSST: %f\n" % (s.Nlosses_T)
            raise exc.NutrientBalanceError(msg)

    def _check_P_balance(self, day):
        s = self.states
        checksum = abs(s.Puptake_T +
                       (self.PamountLVI + self.PamountSTI + self.PamountRTI + self.PamountSOI) -
                       (s.PamountLV + s.PamountST + s.PamountRT + s.PamountSO + s.Plosses_T))

        if abs(checksum) >= 1.:
            msg = "P flows not balanced on day %s\n" % day
            msg += "Checksum: %f, Puptake_T: %f\n" % (checksum, s.Puptake_T)
            msg += "PamountLVI: %f, PamountSTI: %f, PamountRTI: %f, PamountSOI: %f\n" % \
                   (self.PamountLVI, self.PamountSTI, self.PamountRTI, self.PamountSOI)
            msg += "PamountLV: %f, PamountST: %f, PamountRT: %f, PamountSO: %f\n" % \
                   (s.PamountLV, s.PamountST, s.PamountRT, s.PamountSO)
            msg += "PLOSST: %f\n" % (s.Plosses_T)
            raise exc.NutrientBalanceError(msg)

    def _check_K_balance(self, day):
        s = self.states
        checksum = abs(s.Kuptake_T +
                       (self.KamountLVI + self.KamountSTI + self.KamountRTI + self.KamountSOI) -
                       (s.KamountLV + s.KamountST + s.KamountRT + s.KamountSO + s.Klosses_T))

        if abs(checksum) >= 1.:
            msg = "K flows not balanced on day %s\n" % day
            msg += "Checksum: %f, Kuptake_T: %f\n"  % (checksum, s.Kuptake_T)
            msg += "KamountLVI: %f, KamountSTI: %f, KamountRTI: %f, KamountSOI: %f\n" % \
                   (self.KamountLVI, self.KamountSTI, self.KamountRTI, self.KamountSOI)
            msg += "KamountLV: %f, KamountST: %f, KamountRT: %f, KamountSO: %f\n" % \
                   (s.KamountLV, s.KamountST, s.KamountRT, s.KamountSO)
            msg += "KLOSST: %f\n" % (s.Klosses_T)
            raise exc.NutrientBalanceError(msg)
