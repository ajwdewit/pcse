#!/usr/bin/env python

from ...base_classes import StatesTemplate, ParamTemplate, SimulationObject, \
       AfgenTrait
from ...decorators import prepare_rates, prepare_states
from ...traitlets import HasTraits, Float, Int, Instance

class npk_demand(SimulationObject):
    
    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        PMAXLV_TB = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        KMAXLV_TB = AfgenTrait()  # maximum P concentration in leaves as function of dvs
        
        NMAXRT_FR   = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        NMAXST_FR   = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves
        PMAXRT_FR   = Float(-99.)  # maximum P concentration in roots as fraction of maximum P concentration in leaves
        PMAXST_FR   = Float(-99.)  # maximum P concentration in stems as fraction of maximum P concentration in leaves
        KMAXRT_FR   = Float(-99.)  # maximum K concentration in roots as fraction of maximum K concentration in leaves
        KMAXST_FR   = Float(-99.)  # maximum K concentration in stems as fraction of maximum K concentration in leaves
        
        NMAXSO = Float(-99.)  # maximum P concentration in storage organs [kg N kg-1 dry biomass]
        PMAXSO = Float(-99.)  # maximum P concentration in storage organs [kg P kg-1 dry biomass]
        KMAXSO = Float(-99.)  # maximum K concentration in storage organs [kg K kg-1 dry biomass]
        
        TCNT   = Float(-99.)  # time coefficient for N translocation to storage organs [days]
        TCPT   = Float(-99.)  # time coefficient for P translocation to storage organs [days]
        TCKT   = Float(-99.)  # time coefficient for K translocation to storage organs [days]
        
    class StateVariables(StatesTemplate):
        NDEML  = Float(-99.)
        NDEMS  = Float(-99.)
        NDEMR  = Float(-99.)
        NDEMSO = Float(-99.)

        PDEML  = Float(-99.)
        PDEMS  = Float(-99.)
        PDEMR  = Float(-99.)
        PDEMSO = Float(-99.)
        
        KDEML  = Float(-99.)
        KDEMS  = Float(-99.)
        KDEMR  = Float(-99.)
        KDEMSO = Float(-99.)

    
    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk demand __call__()
        """
        self.states = self.StateVariables(kiosk,publish=["NDEML","NDEMS","NDEMR","NDEMSO",
                                                         "PDEML","PDEMS","PDEMR","PDEMSO",
                                                         "KDEML","KDEMS","KDEMR","KDEMSO"],
                                            NDEML=0.,NDEMS=0.,NDEMR=0.,NDEMSO=0.,
                                            PDEML=0.,PDEMS=0.,PDEMR=0.,PDEMSO=0.,
                                            KDEML=0.,KDEMS=0.,KDEMR=0.,KDEMSO=0.)

                                                                                                                
        self.params = self.Parameters(cropdata)
        self.kiosk  = kiosk
    
    @prepare_states
    def integrate(self, day):
        states = self.states
        
        # published states from the kiosk
        DVS = self.kiosk["DVS"]
        WLV = self.kiosk["WLV"]
        WST = self.kiosk["WST"]
        WRT = self.kiosk["WRT"]
        WSO = self.kiosk["WSO"]
        
        ANLV = self.kiosk["ANLV"]
        ANST = self.kiosk["ANST"]
        ANRT = self.kiosk["ANRT"]
        ANSO = self.kiosk["ANSO"]
        
        APLV = self.kiosk["APLV"]
        APST = self.kiosk["APST"]
        APRT = self.kiosk["APRT"]
        APSO = self.kiosk["APSO"]
        
        AKLV = self.kiosk["AKLV"]
        AKST = self.kiosk["AKST"]
        AKRT = self.kiosk["AKRT"]
        AKSO = self.kiosk["AKSO"]
        
        params = self.params

#       Maximum NPK concentrations in leaves [kg N kg-1 DM]
        NMAXLV = params.NMAXLV_TB(DVS)
        PMAXLV = params.PMAXLV_TB(DVS)
        KMAXLV = params.KMAXLV_TB(DVS)
        
#       Maximum NPK concentrations in stems and roots [kg N kg-1 DM]
        NMAXST = params.NMAXST_FR * NMAXLV
        NMAXRT = params.NMAXRT_FR * NMAXLV
        NMAXSO = params.NMAXSO
      
        PMAXST = params.PMAXST_FR * PMAXLV
        PMAXRT = params.PMAXRT_FR * PMAXLV
        PMAXSO = params.PMAXSO      
      
        KMAXST = params.KMAXST_FR * KMAXLV
        KMAXRT = params.KMAXRT_FR * KMAXLV
        KMAXSO = params.KMAXSO

#       N demand [kg ha-1] - maybe should be [kg ha-1 day-1]
        states.NDEML  =  max(NMAXLV*WLV  - ANLV, 0.) # maybe should be divided by one day, see equation 5 Shibu etal 2010
        states.NDEMS  =  max(NMAXST*WST  - ANST, 0.)
        states.NDEMR  =  max(NMAXRT*WRT  - ANRT, 0.)
        states.NDEMSO =  max(NMAXSO*WSO  - ANSO, 0.)/params.TCNT

#       P demand [kg ha-1]
        states.PDEML  =  max(PMAXLV*WLV  - APLV, 0.)
        states.PDEMS  =  max(PMAXST*WST  - APST, 0.)
        states.PDEMR  =  max(PMAXRT*WRT  - APRT, 0.)
        states.PDEMSO =  max(PMAXSO*WSO  - APSO, 0.)/params.TCPT

#       K demand [kg ha-1]
        states.KDEML  =  max(KMAXLV*WLV  - AKLV, 0.)
        states.KDEMS  =  max(KMAXST*WST  - AKST, 0.)
        states.KDEMR  =  max(KMAXRT*WRT  - AKRT, 0.)
        states.KDEMSO =  max(KMAXSO*WSO  - AKSO, 0.)/params.TCKT
        
        