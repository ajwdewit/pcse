#!/usr/bin/env python
"""
Class to calculate various nutrient relates stress factors:
    NNI      nitrogen nutrition index   
    NPKI     NPK nutrition index (=minimum of N/P/K-index)
    NPKREF   assimilation reduction factor  
"""


from ...traitlets import Float, Int, Instance, AfgenTrait
from ...util import limit, doy
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk
from ... import exceptions as exc

class npk_stress(SimulationObject):
    """Implementation of npk stress calculation through [NPK]nutrition index.
    
    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    =======  ============================================= =======  ============

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    =======  ================================================= ==== ============

    
    **External dependencies**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology           -
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):
        NMAXLV = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        PMAXLV = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        KMAXLV = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        FRNX  = Float(-99.)   # optimal N concentration as fraction of maximum N concentration
        FRPX  = Float(-99.)   # optimal P concentration as fraction of maximum P concentration
        FRKX  = Float(-99.)   # optimal K concentration as fraction of maximum K concentration
        LRNR   = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        LSNR   = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves
        LRPR   = Float(-99.)  # maximum P concentration in roots as fraction of maximum P concentration in leaves
        LSPR   = Float(-99.)  # maximum P concentration in stems as fraction of maximum P concentration in leaves
        LRKR   = Float(-99.)  # maximum K concentration in roots as fraction of maximum K concentration in leaves
        LSKR   = Float(-99.)  # maximum K concentration in stems as fraction of maximum K concentration in leaves
        RNFLV  = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        RNFST  = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        RPFLV  = Float(-99.)  # residual P fraction in leaves [kg P kg-1 dry biomass]
        RPFST  = Float(-99.)  # residual P fraction in stems [kg P kg-1 dry biomass]
        RKFLV  = Float(-99.)  # residual K fraction in leaves [kg K kg-1 dry biomass]
        RKFST  = Float(-99.)  # residual K fraction in stems [kg K kg-1 dry biomass]
        NLUE   = Float(-99.)  # coefficient for the reduction of RUE due to nutrient (N-P-K) stress
        

    def initialize(self, day, kiosk, cropdata):
        """
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: nutrient stress using __call__()
        """

        self.params = self.Parameters(cropdata)
        self.kiosk  = kiosk
              
        
    def __call__(self, day):
        params = self.params

        # published states from the kiosk
        WLV  = self.kiosk["WLV"]
        WST  = self.kiosk["WST"]
        DVS  = self.kiosk["DVS"]
        
        ANLV = self.kiosk["ANLV"] # N concentration in leaves [kg ha-1]
        ANST = self.kiosk["ANST"] # N concentration in stems [kg ha-1
        
        APLV = self.kiosk["APLV"] # P concentration in leaves [kg ha-1]
        APST = self.kiosk["APST"] # P concentration in stems [kg ha-1]
        
        AKLV = self.kiosk["AKLV"] # K concentration in leaves [kg ha-1]
        AKST = self.kiosk["AKST"] # K concentration in stems [kg ha-1]
       
#       Maximum NPK concentrations in leaves (kg N kg-1 DM)        
        NMAXLV = params.NMAXLV(DVS)
        PMAXLV = params.PMAXLV(DVS)
        KMAXLV = params.KMAXLV(DVS)

#       Maximum NPK concentrations in stems (kg N kg-1 DM)
        NMAXST = params.LSNR * NMAXLV
        PMAXST = params.LSPR * PMAXLV
        KMAXST = params.LSKR * KMAXLV
        
#       Total vegetative living above-ground biomass (kg DM ha-1)     
        TBGMR = WLV + WST 
      
#       Optimal NPK amount in vegetative above-ground living biomass and its NPK concentration
        NOPTL  = params.FRNX * NMAXLV * WLV
        NOPTS  = params.FRNX * NMAXST * WST
        
        POPTL = params.FRPX * PMAXLV * WLV
        POPTS = params.FRPX * PMAXST * WST

        KOPTL = params.FRKX * KMAXLV * WLV
        KOPTS = params.FRKX * KMAXST * WST
        
#       if above-ground living biomass = 0 then optimum = 0
        if TBGMR > 0.:
            NOPTMR = (NOPTL + NOPTS)/TBGMR
            POPTMR = (POPTL + POPTS)/TBGMR
            KOPTMR = (KOPTL + KOPTS)/TBGMR
        else:
            NOPTMR = POPTMR = KOPTMR = 0.
        
      
#       NPK concentration in total vegetative living per kg above-ground biomass  (kg N/P/K kg-1 DM)
#       if above-ground living biomass = 0 then concentration = 0
        if TBGMR > 0.:
            NFGMR  = (ANLV + ANST)/TBGMR
            PFGMR  = (APLV + APST)/TBGMR
            KFGMR  = (AKLV + AKST)/TBGMR
        else:
            NFGMR = PFGMR = KFGMR = 0.

      
#       Residual NPK concentration in total vegetative living above-ground biomass  (kg N/P/K kg-1 DM)
#       if above-ground living biomass = 0 then residual concentration = 0
        if TBGMR > 0.:
            NRMR = (WLV * params.RNFLV + WST * params.RNFST)/TBGMR
            PRMR = (WLV * params.RPFLV + WST * params.RPFST)/TBGMR
            KRMR = (WLV * params.RKFLV + WST * params.RKFST)/TBGMR
        else:
            NRMR = PRMR = KRMR = 0
            
#              
        if (NOPTMR - NRMR) > 0.:
            NNI = max(0.001,(NFGMR-NRMR)/(NOPTMR-NRMR))
        else:
            NNI = 0.001
            
        if (POPTMR - PRMR) > 0.:
            PNI = max (0.001,(PFGMR-PRMR)/(POPTMR-PRMR))
        else:
           PNI = 0.001
            
        if  (KOPTMR-KRMR) > 0:   
            KNI = max (0.001,(KFGMR-KRMR)/(KOPTMR-KRMR))
        else:
            KNI = 0.001
      
        NPKI = min(NNI,PNI,KNI)    

#       Nutrient reduction factor
        NPKREF= limit(0., 1.0, 1. - params.NLUE*(1.0001-NPKI)**2)
         
        return NNI, NPKI, NPKREF
