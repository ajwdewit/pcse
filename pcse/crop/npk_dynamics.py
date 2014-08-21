#!/usr/bin/env python

from .. import exceptions as exc
from ..traitlets import Float, Int, Instance, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

from .nutrients import npk_uptake_rate as nutrient_uptake_rate
from .nutrients import npk_translocation as nutrient_translocation_rate
from .nutrients import npk_losses as nutrient_losses_rate
from .nutrients import npk_translocatable as nutrient_translocatable
from .nutrients import npk_demand as crop_demand
from .nutrients import npk_supply_storage_organs as supply_storage_organs


class NPK_Crop_Dynamics(SimulationObject):
    """Implementation of npk dynamics.
    
    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    DVSNLT  DVS above which no crop N-P-K uptake occurs
    LRNR    maximum N concentration in roots as fraction
            of maximum N concentration in leaves
    LSNR    maximum N concentration in stems as fraction
            of maximum N concentration in leaves
    LRPR    maximum P concentration in roots as fraction
            of maximum P concentration in leaves
    LSPR    maximum P concentration in stems as fraction
            of maximum P concentration in leaves
    LRKR    maximum K concentration in roots as fraction
            of maximum K concentration in leaves
    LSKR    maximum K concentration in stems as fraction
            of maximum K concentration in leaves
    NMAXLV  maximum N concentration in leaves as function
            of development stage                                    |kg N kg-1 dry biomass|
    PMAXLV  maximum P concentration in leaves as function
            of development stage                                    |kg P kg-1 dry biomass|
    KMAXLV  maximum K concentration in leaves as function
            of development stage                                    |kg K kg-1 dry biomass|
    TDWI    Initial total crop dry weight                  SCr      |kg ha-1|
    SPA     Specific Pod Area                              SCr      |ha kg-1|
    =======  ============================================= =======  ============    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    ANLV     Actual N amount in living leaves                       |kg N ha-1|
    APLV     Actual P amount in living leaves                       |kg P ha-1|
    AKLV     Actual K amount in living leaves                       |kg K ha-1|
        
    ANST     Actual N amount in living stems                        |kg N ha-1|
    APST     Actual P amount in living stems                        |kg P ha-1|
    AKST     Actual K amount in living stems                        |kg K ha-1|

    ANSO     Actual N amount in living storage organs               |kg N ha-1|
    APSO     Actual P amount in living storage organs               |kg P ha-1|
    AKSO     Actual K amount in living storage organs               |kg K ha-1|
    
    ANRT     Actual N amount in living roots                        |kg N ha-1|
    APRT     Actual P amount in living roots                        |kg P ha-1|
    AKRT     Actual K amount in living roots                        |kg K ha-1|
    
    NUPTT    total absorbed N amount                                |kg N ha-1|
    PUPTT    total absorbed P amount                                |kg P ha-1|
    KUPTT    total absorbed K amount                                |kg K ha-1|
    NFIXTT   total biological fixated N amount                      |kg N ha-1|
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    RNLV     Weight increase (N) in leaves                      N   |kg ha-1 d-1|
    RPLV     Weight increase (P) in leaves                      N   |kg ha-1 d-1|
    RKLV     Weight increase (K) in leaves                      N   |kg ha-1 d-1|
    
    RNST     Weight increase (N) in stems                       N   |kg ha-1 d-1|
    RPST     Weight increase (P) in stems                       N   |kg ha-1 d-1|
    RKST     Weight increase (K) in stems                       N   |kg ha-1 d-1|
        
    RNRT     Weight increase (N) in roots                       N   |kg ha-1 d-1|
    RPRT     Weight increase (P) in roots                       N   |kg ha-1 d-1|
    RKRT     Weight increase (K) in roots                       N   |kg ha-1 d-1|
    
    RNSO     Weight increase (N) in storage organs              N   |kg ha-1 d-1|
    RPSO     Weight increase (P) in storage organs              N   |kg ha-1 d-1|
    RKSO     Weight increase (K) in storage organs              N   |kg ha-1 d-1|
           
    NLOSSR   N loss due to senescence                               |kg ha-1 d-1| 
    PLOSSR   P loss due to senescence                               |kg ha-1 d-1|
    KLOSSR   K loss due to senescence                               |kg ha-1 d-1|

    
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology           -
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    =======  =================================== =================  ============
    """

#   rate calculation components
    uptake_rate = Instance(SimulationObject)
    transl_rate = Instance(SimulationObject)
    dying_rate = Instance(SimulationObject)

#   state calculation components    
    transl_state = Instance(SimulationObject)
    demand = Instance(SimulationObject)
    supply = Instance(SimulationObject)
    
    ANLVI = Float(-99.) # initial soil N amount in leaves 
    ANSTI = Float(-99.) # initial soil N amount in stems
    ANRTI = Float(-99.) # initial soil N amount in roots
    ANSOI = Float(-99.) # initial soil N amount in storage organs
    
    APLVI = Float(-99.) # initial soil P amount in leaves 
    APSTI = Float(-99.) # initial soil P amount in stems
    APRTI = Float(-99.) # initial soil P amount in roots
    APSOI = Float(-99.) # initial soil P amount in storage organs

    AKLVI = Float(-99.) # initial soil K amount in leaves 
    AKSTI = Float(-99.) # initial soil K amount in stems
    AKRTI = Float(-99.) # initial soil K amount in roots
    AKSOI = Float(-99.) # initial soil K amount in storage organs


    class Parameters(ParamTemplate):
        NMAXLV = AfgenTrait()
        PMAXLV = AfgenTrait()
        KMAXLV = AfgenTrait()
        TDWI = Float(-99.)
        LSNR = Float(-99.)
        LRNR = Float(-99.)
        LSPR = Float(-99.)
        LRPR = Float(-99.)
        LSKR = Float(-99.)
        LRKR = Float(-99.)
        DVSNLT = Float(-99.)

    class StateVariables(StatesTemplate):
        ANLV = Float(-99.) # N amount in leaves [kg N ha-1]
        APLV = Float(-99.) # P amount in leaves [kg P ]
        AKLV = Float(-99.) # K amount in leaves [kg K ]
        
        ANST = Float(-99.) # N amount in stems [kg N ]
        APST = Float(-99.) # P amount in stems [kg P ]
        AKST = Float(-99.) # K amount in stems [kg K ]
      
        ANSO = Float(-99.) # N amount in storage organs [kg N ]
        APSO = Float(-99.) # P amount in storage organs [kg P ]
        AKSO = Float(-99.) # K amount in storage organs [kg K ]
        
        ANRT = Float(-99.) # N amount in roots [kg N ]
        APRT = Float(-99.) # P amount in roots [kg P ]
        AKRT = Float(-99.) # K amount in roots [kg K ]
        
        NUPTT = Float(-99.) # total absorbed N amount [kg N ]
        PUPTT = Float(-99.) # total absorbed P amount [kg P ]
        KUPTT = Float(-99.) # total absorbed K amount [kg K ]
        NFIXTT = Float(-99.) # total biological fixated N amount [kg N ]
        
        NLOSST = Float(-99.)
        PLOSST = Float(-99.)
        KLOSST = Float(-99.)

    class RateVariables(RatesTemplate):
        RNLV = Float(-99.)
        RPLV = Float(-99.)
        RKLV = Float(-99.)
        
        RNST = Float(-99.)
        RPST = Float(-99.)
        RKST = Float(-99.)
               
        RNRT = Float(-99.)
        RPRT = Float(-99.)
        RKRT = Float(-99.)
        
        RNSO = Float(-99.)
        RPSO = Float(-99.)
        RKSO = Float(-99.)
               
        NLOSSR = Float(-99.)
        PLOSSR = Float(-99.)
        KLOSSR = Float(-99.)
        
    def initialize(self, day, kiosk, cropdata):
        """
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """  
        
        self.params = self.Parameters(cropdata)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
#       Initialize components of the npk_crop_dynamics
        self.uptake_rate = nutrient_uptake_rate(day, kiosk, cropdata)
        self.transl_rate = nutrient_translocation_rate(day, kiosk, cropdata)
        self.dying_rate  = nutrient_losses_rate(day, kiosk, cropdata)
        
        self.transl_state = nutrient_translocatable(day, kiosk, cropdata)
        self.demand = crop_demand(day, kiosk, cropdata)
        self.supply = supply_storage_organs(day, kiosk, cropdata)
        
        # INITIAL STATES
        params = self.params
        states = self.states
        # Initial storage organ biomass
        FO   = self.kiosk["FO"]
        FR   = self.kiosk["FR"]
        FS   = self.kiosk["FS"]
        FL   = self.kiosk["FL"]
        
        DVS  = self.kiosk["DVS"]
               
        WLV  = params.TDWI * (1-FR) * FL
        WST  = params.TDWI * (1-FR) * FS
        WSO  = params.TDWI * (1-FR) * FO
        WRT  = params.TDWI * FR
        
        self.ANLVI = ANLV = WLV * params.NMAXLV(DVS)
        self.ANSTI = ANST = WST * params.NMAXLV(DVS) * params.LSNR
        self.ANRTI = ANRT = WRT * params.NMAXLV(DVS) * params.LRNR
        self.ANSOI = ANSO = 0.
        
        self.APLVI = APLV = WLV * params.PMAXLV(DVS)
        self.APSTI = APST = WST * params.PMAXLV(DVS) * params.LSPR
        self.APRTI = APRT = WRT * params.PMAXLV(DVS) * params.LRPR
        self.APSOI = APSO = 0.

        self.AKLVI = AKLV = WLV * params.KMAXLV(DVS)
        self.AKSTI = AKST = WST * params.KMAXLV(DVS) * params.LSKR
        self.AKRTI = AKRT = WRT * params.KMAXLV(DVS) * params.LRKR
        self.AKSOI = AKSO = 0.
        
        params = self.params
        DVSNLT = params.DVSNLT

        self.states = self.StateVariables(kiosk,
                        publish=["ANLV","ANST","ANRT","ANSO", "APLV","APST",
                                 "APRT","APSO", "AKLV","AKST","AKRT","AKSO"],
                        ANLV=ANLV, ANST=ANST, ANRT=ANRT, ANSO=ANSO,
                        APLV=APLV, APST=APST, APRT=APRT, APSO=APSO,
                        AKLV=AKLV, AKST=AKST, AKRT=AKRT, AKSO=AKSO,
                        NUPTT=0 ,PUPTT=0., KUPTT=0., NFIXTT=0.,
                        NLOSST=0  ,PLOSST=0., KLOSST=0.)
        
   
    @staticmethod
    def _check_N_balance(self, day):
        states = self.states
        
        NUPTT  = states.NUPTT
        NFIXTT = states.NFIXTT
        
        ANLVI  = self.ANLVI
        ANSTI  = self.ANSTI
        ANRTI  = self.ANRTI
        ANSOI  = self.ANSOI
        
        ANLV  = states.ANLV
        ANST  = states.ANST
        ANRT  = states.ANRT
        ANSO  = states.ANSO
        
        NLOSST = states.NLOSST

        checksum = abs(NUPTT + NFIXTT + (ANLVI + ANSTI + ANRTI + ANSOI) - \
                    (ANLV + ANST + ANRT + ANSO + NLOSST))
        
        if abs(checksum) >= 1.:
            msg = "N flows not balanced on day %s\n" % day
            msg += "Checksum: %f, NUPTT: %f, NFIXTT: %f\n" % (checksum, NUPTT, NFIXTT)
            msg += "ANLVI: %f, ANSTI: %f, ANRTI: %f, ANSOI: %f\n"  %(ANLVI, ANSTI, ANRTI, ANSOI)
            msg += "ANLV: %f, ANST: %f, ANRT: %f, ANSO: %f\n" % (ANLV, ANST, ANRT, ANSO) 
            msg += "NLOSST: %f\n" %(NLOSST)
            raise exc.NutrientBalanceError(msg)
            
     
    @staticmethod
    def _check_P_balance(self, day):
        states = self.states            
        PUPTT  = states.PUPTT
             
        APLVI  = self.APLVI
        APSTI  = self.APSTI
        APRTI  = self.APRTI
        APSOI  = self.APSOI
        
        APLV  = states.APLV
        APST  = states.APST
        APRT  = states.APRT
        APSO  = states.APSO
        
        PLOSST = states.PLOSST
        
        checksum = abs(PUPTT + (APLVI + APSTI + APRTI + APSOI) - \
                    (APLV + APST + APRT + APSO + PLOSST))
            
        if abs(checksum) >= 1.:
            msg = "P flows not balanced on day %s\n" % day
            msg += "Checksum: %f, PUPTT: %f\n" % (checksum, PUPTT)
            msg += "APLVI: %f, APSTI: %f, APRTI: %f, APSOI: %f\n" % (APLVI, APSTI, APRTI, APSOI)
            msg += "APLV: %f, APST: %f, APRT: %f, APSO: %f\n" % (APLV, APST, APRT, APSO) 
            msg += "PLOSST: %f\n" %(PLOSST)
            raise exc.NutrientBalanceError(msg)
            
    @staticmethod
    def _check_K_balance(self, day):
        states = self.states            
        KUPTT  = states.KUPTT
             
        AKLVI  = self.AKLVI
        AKSTI  = self.AKSTI
        AKRTI  = self.AKRTI
        AKSOI  = self.AKSOI
        
        AKLV  = states.AKLV
        AKST  = states.AKST
        AKRT  = states.AKRT
        AKSO  = states.AKSO
        
        KLOSST = states.KLOSST
        
        checksum = abs(KUPTT + (AKLVI + AKSTI + AKRTI + AKSOI) - \
                    (AKLV + AKST + AKRT + AKSO + KLOSST))
            
        if abs(checksum) >= 1.:
            msg = "K flows not balanced on day %s\n" % day
            msg += "Checksum: %f, KUPTT: %f\n"  %(checksum, KUPTT)
            msg += "AKLVI: %f, AKSTI: %f, AKRTI: %f, AKSOI: %f\n" % (AKLVI, AKSTI, AKRTI, AKSOI)
            msg += "AKLV: %f, AKST: %f, AKRT: %f, AKSO: %f\n" % (AKLV, AKST, AKRT, AKSO) 
            msg += "KLOSST: %f\n" %(KLOSST)
            raise exc.NutrientBalanceError(msg)    
            

    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates     
        
        self.uptake_rate.calc_rates(day)
        self.transl_rate.calc_rates(day)
        self.dying_rate.calc_rates(day)
        
                   
        # N rates in leaves, stems, root and storage organs
        rates.RNLV = self.kiosk["RNULV"] - self.kiosk["RNTLV"] - self.kiosk["RNDLV"]
        rates.RNST = self.kiosk["RNUST"] - self.kiosk["RNTST"] - self.kiosk["RNDST"]
        rates.RNRT = self.kiosk["RNURT"] - self.kiosk["RNTRT"] - self.kiosk["RNDRT"]
        rates.RNSO = self.kiosk["RNUSO"]
        
        # P rates in leaves, stems, root and storage organs
        rates.RPLV = self.kiosk["RPULV"] - self.kiosk["RPTLV"] - self.kiosk["RPDLV"]
        rates.RPST = self.kiosk["RPUST"] - self.kiosk["RPTST"] - self.kiosk["RPDST"]
        rates.RPRT = self.kiosk["RPURT"] - self.kiosk["RPTRT"] - self.kiosk["RPDRT"]
        rates.RPSO = self.kiosk["RPUSO"]

        # K rates in leaves, stems, root and storage organs
        rates.RKLV = self.kiosk["RKULV"] - self.kiosk["RKTLV"] - self.kiosk["RKDLV"]
        rates.RKST = self.kiosk["RKUST"] - self.kiosk["RKTST"] - self.kiosk["RKDST"]
        rates.RKRT = self.kiosk["RKURT"] - self.kiosk["RKTRT"] - self.kiosk["RKDRT"]
        rates.RKSO = self.kiosk["RKUSO"]
        
        rates.NLOSSR = self.kiosk["RNDLV"] + self.kiosk["RNDST"] + self.kiosk["RNDRT"]
        rates.PLOSSR = self.kiosk["RPDLV"] + self.kiosk["RPDST"] + self.kiosk["RPDRT"]
        rates.KLOSSR = self.kiosk["RKDLV"] + self.kiosk["RKDST"] + self.kiosk["RKDRT"] 

        self._check_N_balance(self, day)
        self._check_P_balance(self, day)
        self._check_K_balance(self, day)
        
        
    @prepare_states
    def integrate(self, day):
        rates  = self.rates
        states = self.states
        

        # N amount in leaves, stems, root and storage organs
        states.ANLV += rates.RNLV
        states.ANST += rates.RNST
        states.ANRT += rates.RNRT
        states.ANSO += rates.RNSO
        
        # P amount in leaves, stems, root and storage organs
        states.APLV += rates.RPLV
        states.APST += rates.RPST
        states.APRT += rates.RPRT
        states.APSO += rates.RPSO

        # K amount in leaves, stems, root and storage organs
        states.AKLV += rates.RKLV
        states.AKST += rates.RKST
        states.AKRT += rates.RKRT
        states.AKSO += rates.RKSO
        
        # translocatable NPK amount
        self.transl_state.integrate(day)
        # NPK demand
        self.demand.integrate(day)
        # NPK supply to storage organs
        self.supply.integrate(day)
        
        # total NPK uptake from soil
        states.NUPTT  += self.kiosk["NUPTR"]
        states.PUPTT  += self.kiosk["PUPTR"]
        states.KUPTT  += self.kiosk["KUPTR"]
        states.NFIXTT += self.kiosk["NFIXTR"]
        
        states.NLOSST += rates.NLOSSR
        states.PLOSST += rates.PLOSSR
        states.KLOSST += rates.KLOSSR
    
    
    