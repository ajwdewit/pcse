#!/usr/bin/env python

from .. import exceptions as exc
from ..traitlets import Float, Int, Instance, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

#from .nutrients import npk_uptake_rate as nutrient_uptake_rate
#from .nutrients import npk_translocation as nutrient_translocation_rate
from .nutrients import npk_losses as nutrient_losses_rate
from .nutrients import npk_translocatable as nutrient_translocatable
from .nutrients import NPK_Demand_Uptake as crop_demand
#from .nutrients import npk_supply_storage_organs as supply_storage_organs


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
    
    NUPTAKE_T    total absorbed N amount                                |kg N ha-1|
    PUPTAKE_T    total absorbed P amount                                |kg P ha-1|
    KUPTAKE_T    total absorbed K amount                                |kg K ha-1|
    NFIX_T   total biological fixated N amount                      |kg N ha-1|
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
#    uptake_rate = Instance(SimulationObject)
    translocation = Instance(SimulationObject)
    dying_rate = Instance(SimulationObject)

#   state calculation components    
#    transl_state = Instance(SimulationObject)
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
        NMAXLV_TB = AfgenTrait()
        PMAXLV_TB = AfgenTrait()
        KMAXLV_TB = AfgenTrait()
        TDWI = Float(-99.)
        NMAXST_FR = Float(-99.)
        NMAXRT_FR = Float(-99.)
        PMAXST_FR = Float(-99.)
        PMAXRT_FR = Float(-99.)
        KMAXST_FR = Float(-99.)
        KMAXRT_FR = Float(-99.)
        DVSNPK_STOP = Float(-99.)

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
        
        NUPTAKE_T = Float(-99.) # total absorbed N amount [kg N ]
        PUPTAKE_T = Float(-99.) # total absorbed P amount [kg P ]
        KUPTAKE_T = Float(-99.) # total absorbed K amount [kg K ]
        NFIX_T = Float(-99.) # total biological fixated N amount [kg N ]
        
        NLOSSES_T = Float(-99.)
        PLOSSES_T = Float(-99.)
        KLOSSES_T = Float(-99.)

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
               
        RNLOSS = Float(-99.)
        RPLOSS = Float(-99.)
        RKLOSS = Float(-99.)
        
    def initialize(self, day, kiosk, cropdata):
        """
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """  
        
        self.params = self.Parameters(cropdata)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
#       Initialize components of the npk_crop_dynamics
#        self.uptake_rate = nutrient_uptake_rate(day, kiosk, cropdata)
#        self.transl_rate = nutrient_translocation_rate(day, kiosk, cropdata)
        self.dying_rate  = nutrient_losses_rate(day, kiosk, cropdata)
        
        self.translocation = nutrient_translocatable(day, kiosk, cropdata)
        self.demand = crop_demand(day, kiosk, cropdata)
#        self.supply = supply_storage_organs(day, kiosk, cropdata)
        
        # INITIAL STATES
        params = self.params
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
        
        self.ANLVI = ANLV = WLV * params.NMAXLV_TB(DVS)
        self.ANSTI = ANST = WST * params.NMAXLV_TB(DVS) * params.NMAXST_FR
        self.ANRTI = ANRT = WRT * params.NMAXLV_TB(DVS) * params.NMAXRT_FR
        self.ANSOI = ANSO = 0.
        
        self.APLVI = APLV = WLV * params.PMAXLV_TB(DVS)
        self.APSTI = APST = WST * params.PMAXLV_TB(DVS) * params.PMAXST_FR
        self.APRTI = APRT = WRT * params.PMAXLV_TB(DVS) * params.PMAXRT_FR
        self.APSOI = APSO = 0.

        self.AKLVI = AKLV = WLV * params.KMAXLV_TB(DVS)
        self.AKSTI = AKST = WST * params.KMAXLV_TB(DVS) * params.KMAXST_FR
        self.AKRTI = AKRT = WRT * params.KMAXLV_TB(DVS) * params.KMAXRT_FR
        self.AKSOI = AKSO = 0.

        self.states = self.StateVariables(kiosk,
                        publish=["ANLV","ANST","ANRT","ANSO", "APLV","APST",
                                 "APRT","APSO", "AKLV","AKST","AKRT","AKSO"],
                        ANLV=ANLV, ANST=ANST, ANRT=ANRT, ANSO=ANSO,
                        APLV=APLV, APST=APST, APRT=APRT, APSO=APSO,
                        AKLV=AKLV, AKST=AKST, AKRT=AKRT, AKSO=AKSO,
                        NUPTAKE_T=0 ,PUPTAKE_T=0., KUPTAKE_T=0., NFIX_T=0.,
                        NLOSSES_T=0 ,PLOSSES_T=0., KLOSSES_T=0.)


    @prepare_rates
    def calc_rates(self, day):
        rates = self.rates
        
        self.demand.calc_rates(day)
        self.translocation.calc_rates(day)
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
        
        rates.RNLOSS = self.kiosk["RNDLV"] + self.kiosk["RNDST"] + self.kiosk["RNDRT"]
        rates.RPLOSS = self.kiosk["RPDLV"] + self.kiosk["RPDST"] + self.kiosk["RPDRT"]
        rates.RKLOSS = self.kiosk["RKDLV"] + self.kiosk["RKDST"] + self.kiosk["RKDRT"]

        self._check_N_balance(day)
        self._check_P_balance(day)
        self._check_K_balance(day)
        
        
    @prepare_states
    def integrate(self, day):
        rates = self.rates
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
        self.translocation.integrate(day)
        # NPK demand
        self.demand.integrate(day)
        # NPK supply to storage organs
#        self.supply.integrate(day)
        
        # total NPK uptake from soil
        states.NUPTAKE_T += self.kiosk["RNUPTAKE"]
        states.PUPTAKE_T += self.kiosk["RPUPTAKE"]
        states.KUPTAKE_T += self.kiosk["RKUPTAKE"]
        states.NFIX_T += self.kiosk["RNFIX"]
        
        states.NLOSSES_T += rates.RNLOSS
        states.PLOSSES_T += rates.RPLOSS
        states.KLOSSES_T += rates.RKLOSS

    def _check_N_balance(self, day):
        states = self.states

        NUPTAKE_T  = states.NUPTAKE_T
        NFIX_T = states.NFIX_T

        ANLVI  = self.ANLVI
        ANSTI  = self.ANSTI
        ANRTI  = self.ANRTI
        ANSOI  = self.ANSOI

        ANLV  = states.ANLV
        ANST  = states.ANST
        ANRT  = states.ANRT
        ANSO  = states.ANSO

        NLOSST = states.NLOSSES_T

        checksum = abs(NUPTAKE_T + NFIX_T + (ANLVI + ANSTI + ANRTI + ANSOI) -
                       (ANLV + ANST + ANRT + ANSO + NLOSST))

        if abs(checksum) >= 1.:
            msg = "N flows not balanced on day %s\n" % day
            msg += "Checksum: %f, NUPTAKE_T: %f, NFIX_T: %f\n" % (checksum, NUPTAKE_T, NFIX_T)
            msg += "ANLVI: %f, ANSTI: %f, ANRTI: %f, ANSOI: %f\n"  %(ANLVI, ANSTI, ANRTI, ANSOI)
            msg += "ANLV: %f, ANST: %f, ANRT: %f, ANSO: %f\n" % (ANLV, ANST, ANRT, ANSO)
            msg += "NLOSST: %f\n" %(NLOSST)
            raise exc.NutrientBalanceError(msg)


    def _check_P_balance(self, day):
        states = self.states
        PUPTAKE_T = states.PUPTAKE_T

        APLVI  = self.APLVI
        APSTI  = self.APSTI
        APRTI  = self.APRTI
        APSOI  = self.APSOI

        APLV  = states.APLV
        APST  = states.APST
        APRT  = states.APRT
        APSO  = states.APSO

        PLOSST = states.PLOSSES_T

        checksum = abs(PUPTAKE_T + (APLVI + APSTI + APRTI + APSOI) - \
                    (APLV + APST + APRT + APSO + PLOSST))

        if abs(checksum) >= 1.:
            msg = "P flows not balanced on day %s\n" % day
            msg += "Checksum: %f, PUPTAKE_T: %f\n" % (checksum, PUPTAKE_T)
            msg += "APLVI: %f, APSTI: %f, APRTI: %f, APSOI: %f\n" % (APLVI, APSTI, APRTI, APSOI)
            msg += "APLV: %f, APST: %f, APRT: %f, APSO: %f\n" % (APLV, APST, APRT, APSO)
            msg += "PLOSST: %f\n" %(PLOSST)
            raise exc.NutrientBalanceError(msg)

    def _check_K_balance(self, day):
        states = self.states
        KUPTAKE_T  = states.KUPTAKE_T

        AKLVI  = self.AKLVI
        AKSTI  = self.AKSTI
        AKRTI  = self.AKRTI
        AKSOI  = self.AKSOI

        AKLV  = states.AKLV
        AKST  = states.AKST
        AKRT  = states.AKRT
        AKSO  = states.AKSO

        KLOSST = states.KLOSSES_T

        checksum = abs(KUPTAKE_T + (AKLVI + AKSTI + AKRTI + AKSOI) - \
                    (AKLV + AKST + AKRT + AKSO + KLOSST))

        if abs(checksum) >= 1.:
            msg = "K flows not balanced on day %s\n" % day
            msg += "Checksum: %f, KUPTAKE_T: %f\n"  %(checksum, KUPTAKE_T)
            msg += "AKLVI: %f, AKSTI: %f, AKRTI: %f, AKSOI: %f\n" % (AKLVI, AKSTI, AKRTI, AKSOI)
            msg += "AKLV: %f, AKST: %f, AKRT: %f, AKSO: %f\n" % (AKLV, AKST, AKRT, AKSO)
            msg += "KLOSST: %f\n" %(KLOSST)
            raise exc.NutrientBalanceError(msg)
