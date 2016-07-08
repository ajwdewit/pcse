#!/usr/bin/env python

from .. import exceptions as exc
from ..traitlets import Float, Int, Instance, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

from .nutrients import NPK_Translocation
from .nutrients import NPK_Demand_Uptake

class NPK_Crop_Dynamics(SimulationObject):
    """Implementation of overall NPK crop dynamics.

    NPK_Crop_Dynamics implements the overall logic of N/P/K book-keeping within the
    crop.

    **Simulation parameters**
    
    ===========  ================================================ =======  ======================
     Name         Description                                      Type     Unit
    ===========  ================================================ =======  ======================
    DVSNPK_STOP   DVS above which no crop N-P-K uptake occurs      SCr       -

    NMAXLV_TB      Maximum N concentration in leaves as            TCr     kg N kg-1 dry biomass
                   function of dvs
    PMAXLV_TB      As for P                                        TCr     kg P kg-1 dry biomass
    KMAXLV_TB      As for K                                        TCr     kg K kg-1 dry biomass

    NMAXRT_FR      Maximum N concentration in roots as fraction    SCr     -
                   of maximum N concentration in leaves
    PMAXRT_FR      As for P                                        SCr     -
    KMAXRT_FR      As for K                                        SCr     -

    NMAXST_FR      Maximum N concentration in stems as fraction    SCr     -
                   of maximum N concentration in leaves
    KMAXST_FR      As for K                                        SCr     -
    PMAXST_FR      As for P                                        SCr     -

    NRESIDLV       Residual N fraction in leaves                   SCr     kg N kg-1 dry biomass
    PRESIDLV       Residual P fraction in leaves                   SCr     kg P kg-1 dry biomass
    KRESIDLV       Residual K fraction in leaves                   SCr     kg K kg-1 dry biomass

    NRESIDRT       Residual N fraction in roots                    SCr     kg N kg-1 dry biomass
    PRESIDRT       Residual P fraction in roots                    SCr     kg P kg-1 dry biomass
    KRESIDRT       Residual K fraction in roots                    SCr     kg K kg-1 dry biomass

    NRESIDST       Residual N fraction in stems                    SCr     kg N kg-1 dry biomass
    PRESIDST       Residual P fraction in stems                    SCr     kg P kg-1 dry biomass
    KRESIDST       Residual K fraction in stems                    SCr     kg K kg-1 dry biomass
    ===========  ================================================ =======  ======================

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    ANLV     Actual N amount in living leaves                   Y   |kg N ha-1|
    APLV     Actual P amount in living leaves                   Y   |kg P ha-1|
    AKLV     Actual K amount in living leaves                   Y   |kg K ha-1|
        
    ANST     Actual N amount in living stems                    Y   |kg N ha-1|
    APST     Actual P amount in living stems                    Y   |kg P ha-1|
    AKST     Actual K amount in living stems                    Y   |kg K ha-1|

    ANSO     Actual N amount in living storage organs           Y   |kg N ha-1|
    APSO     Actual P amount in living storage organs           Y   |kg P ha-1|
    AKSO     Actual K amount in living storage organs           Y   |kg K ha-1|
    
    ANRT     Actual N amount in living roots                    Y   |kg N ha-1|
    APRT     Actual P amount in living roots                    Y   |kg P ha-1|
    AKRT     Actual K amount in living roots                    Y   |kg K ha-1|
    
    NUPTAKE_T    total absorbed N amount                        N   |kg N ha-1|
    PUPTAKE_T    total absorbed P amount                        N   |kg P ha-1|
    KUPTAKE_T    total absorbed K amount                        N   |kg K ha-1|
    NFIX_T   total biological fixated N amount                  N   |kg N ha-1|
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

    RNDLV    Rate of N loss in leaves                           N   |kg ha-1 d-1|
    RPDLV    as for P                                           N   |kg ha-1 d-1|
    RKDLV    as for K                                           N   |kg ha-1 d-1|

    RNDST    Rate of N loss in roots                            N   |kg ha-1 d-1|
    RPDST    as for P                                           N   |kg ha-1 d-1|
    RKDST    as for K                                           N   |kg ha-1 d-1|

    RNDRT    Rate of N loss in stems                            N   |kg ha-1 d-1|
    RPDRT    as for P                                           N   |kg ha-1 d-1|
    RKDRT    as for K                                           N   |kg ha-1 d-1|

    RNLOSS   N loss due to senescence                           N   |kg ha-1 d-1|
    RPLOSS   P loss due to senescence                           N   |kg ha-1 d-1|
    RKLOSS   K loss due to senescence                           N   |kg ha-1 d-1|
    =======  ================================================= ==== ============
    
    **Signals send or handled**
    
    None
    
    **External dependencies**
    
    =======  =================================== ====================  ============
     Name     Description                         Provided by            Unit
    =======  =================================== ====================  ============
    DVS      Crop development stage              DVS_Phenology           -
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    WRT      Dry weight of living roots          WOFOST_Root_Dynamics  |kg ha-1|
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    DRLV     Death rate of leaves                WOFOST_Leaf_Dynamics  |kg ha-1 d-|
    DRRT     Death rate of roots                 WOFOST_Root_Dynamics  |kg ha-1 d-|
    DRST     Death rate of stems                 WOFOST_Stem_Dynamics  |kg ha-1 d-|
    =======  =================================== ====================  ============
    """

    translocation = Instance(SimulationObject)
    demand_uptake = Instance(SimulationObject)

    ANLVI = Float(-99.)  # initial soil N amount in leaves
    ANSTI = Float(-99.)  # initial soil N amount in stems
    ANRTI = Float(-99.)  # initial soil N amount in roots
    ANSOI = Float(-99.)  # initial soil N amount in storage organs
    
    APLVI = Float(-99.)  # initial soil P amount in leaves
    APSTI = Float(-99.)  # initial soil P amount in stems
    APRTI = Float(-99.)  # initial soil P amount in roots
    APSOI = Float(-99.)  # initial soil P amount in storage organs

    AKLVI = Float(-99.)  # initial soil K amount in leaves
    AKSTI = Float(-99.)  # initial soil K amount in stems
    AKRTI = Float(-99.)  # initial soil K amount in roots
    AKSOI = Float(-99.)  # initial soil K amount in storage organs

    class Parameters(ParamTemplate):
        DVSNPK_STOP = Float(-99.)
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
        RNLV = Float(-99.)  # Net rates of NPK in different plant organs 
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
               
        RNDLV = Float(-99.)  # N loss rate leaves [kg ha-1 d-1]
        RNDST = Float(-99.)  # N loss rate stems  [kg ha-1 d-1]
        RNDRT = Float(-99.)  # N loss rate roots  [kg ha-1 d-1]
        
        RPDLV = Float(-99.)  # P loss rate leaves [kg ha-1 d-1]
        RPDST = Float(-99.)  # P loss rate stems  [kg ha-1 d-1]
        RPDRT = Float(-99.)  # P loss rate roots  [kg ha-1 d-1]
        
        RKDLV = Float(-99.)  # K loss rate leaves [kg ha-1 d-1]
        RKDST = Float(-99.)  # K loss rate stems  [kg ha-1 d-1]
        RKDRT = Float(-99.)  # K loss rate roots  [kg ha-1 d-1]

        RNLOSS = Float(-99.)
        RPLOSS = Float(-99.)
        RKLOSS = Float(-99.)
        
    def initialize(self, day, kiosk, parvalues):
        """
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameters as key/value pairs
        """  
        
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
#       Initialize components of the npk_crop_dynamics
        self.translocation = NPK_Translocation(day, kiosk, parvalues)
        self.demand_uptake = NPK_Demand_Uptake(day, kiosk, parvalues)

        # INITIAL STATES
        params = self.params

        DVS = self.kiosk["DVS"]
        WLV = self.kiosk["WLV"]
        WST = self.kiosk["WST"]
        WRT = self.kiosk["WRT"]

        # Initial amounts
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
                        publish=["ANLV", "ANST", "ANRT", "ANSO", "APLV", "APST",
                                 "APRT", "APSO", "AKLV", "AKST", "AKRT", "AKSO"],
                        ANLV=ANLV, ANST=ANST, ANRT=ANRT, ANSO=ANSO,
                        APLV=APLV, APST=APST, APRT=APRT, APSO=APSO,
                        AKLV=AKLV, AKST=AKST, AKRT=AKRT, AKSO=AKSO,
                        NUPTAKE_T=0, PUPTAKE_T=0., KUPTAKE_T=0., NFIX_T=0.,
                        NLOSSES_T=0, PLOSSES_T=0., KLOSSES_T=0.)

    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        params = self.params
        
        self.demand_uptake.calc_rates(day, drv)
        self.translocation.calc_rates(day, drv)

        # Compute loss of NPK due to death of plant material
        DRLV = self.kiosk["DRLV"]  # death rate leaves [kg dry matter ha-1 d-1]
        DRST = self.kiosk["DRST"]  # death rate stems [kg dry matter ha-1 d-1]
        DRRT = self.kiosk["DRRT"]  # death rate roots [kg dry matter ha-1 d-1]

        rates.RNDLV = params.NRESIDLV * DRLV
        rates.RNDST = params.NRESIDST * DRST
        rates.RNDRT = params.NRESIDRT * DRRT

        rates.RPDLV = params.PRESIDLV * DRLV
        rates.RPDST = params.PRESIDST * DRST
        rates.RPDRT = params.PRESIDRT * DRRT

        rates.RKDLV = params.KRESIDLV * DRLV
        rates.RKDST = params.KRESIDST * DRST
        rates.RKDRT = params.KRESIDRT * DRRT

        # N rates in leaves, stems, root and storage organs computed as
        # uptake - translocation - death.
        # except for storage organs which only take up as a result of translocation.
        rates.RNLV = self.kiosk["RNULV"] - self.kiosk["RNTLV"] - rates.RNDLV
        rates.RNST = self.kiosk["RNUST"] - self.kiosk["RNTST"] - rates.RNDST
        rates.RNRT = self.kiosk["RNURT"] - self.kiosk["RNTRT"] - rates.RNDRT
        rates.RNSO = self.kiosk["RNUSO"]
        
        # P rates in leaves, stems, root and storage organs
        rates.RPLV = self.kiosk["RPULV"] - self.kiosk["RPTLV"] - rates.RPDLV
        rates.RPST = self.kiosk["RPUST"] - self.kiosk["RPTST"] - rates.RPDST
        rates.RPRT = self.kiosk["RPURT"] - self.kiosk["RPTRT"] - rates.RPDRT
        rates.RPSO = self.kiosk["RPUSO"]

        # K rates in leaves, stems, root and storage organs
        rates.RKLV = self.kiosk["RKULV"] - self.kiosk["RKTLV"] - rates.RKDLV
        rates.RKST = self.kiosk["RKUST"] - self.kiosk["RKTST"] - rates.RKDST
        rates.RKRT = self.kiosk["RKURT"] - self.kiosk["RKTRT"] - rates.RKDRT
        rates.RKSO = self.kiosk["RKUSO"]
        
        rates.RNLOSS = rates.RNDLV + rates.RNDST + rates.RNDRT
        rates.RPLOSS = rates.RPDLV + rates.RPDST + rates.RPDRT
        rates.RKLOSS = rates.RKDLV + rates.RKDST + rates.RKDRT

        self._check_N_balance(day)
        self._check_P_balance(day)
        self._check_K_balance(day)
        
    @prepare_states
    def integrate(self, day, delt=1.0):
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
        
        self.translocation.integrate(day, delt)
        self.demand_uptake.integrate(day, delt)

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

        APLVI = self.APLVI
        APSTI = self.APSTI
        APRTI = self.APRTI
        APSOI = self.APSOI

        APLV = states.APLV
        APST = states.APST
        APRT = states.APRT
        APSO = states.APSO

        PLOSST = states.PLOSSES_T

        checksum = abs(PUPTAKE_T + (APLVI + APSTI + APRTI + APSOI) -
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
        KUPTAKE_T = states.KUPTAKE_T

        AKLVI = self.AKLVI
        AKSTI = self.AKSTI
        AKRTI = self.AKRTI
        AKSOI = self.AKSOI

        AKLV = states.AKLV
        AKST = states.AKST
        AKRT = states.AKRT
        AKSO = states.AKSO

        KLOSST = states.KLOSSES_T

        checksum = abs(KUPTAKE_T + (AKLVI + AKSTI + AKRTI + AKSOI) -
                       (AKLV + AKST + AKRT + AKSO + KLOSST))

        if abs(checksum) >= 1.:
            msg = "K flows not balanced on day %s\n" % day
            msg += "Checksum: %f, KUPTAKE_T: %f\n"  %(checksum, KUPTAKE_T)
            msg += "AKLVI: %f, AKSTI: %f, AKRTI: %f, AKSOI: %f\n" % (AKLVI, AKSTI, AKRTI, AKSOI)
            msg += "AKLV: %f, AKST: %f, AKRT: %f, AKSO: %f\n" % (AKLV, AKST, AKRT, AKSO)
            msg += "KLOSST: %f\n" %(KLOSST)
            raise exc.NutrientBalanceError(msg)
