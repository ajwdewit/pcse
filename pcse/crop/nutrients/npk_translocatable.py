#!/usr/bin/env python

from ...traitlets import Float, Instance
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class npk_translocatable(SimulationObject):
    
    class Parameters(ParamTemplate):
        RNFLV = Float(-99.) # residual N fraction in leaves [kg N kg-1 dry biomass]
        RNFST = Float(-99.) # residual N fraction in stems [kg N kg-1 dry biomass]
        RNFRT = Float(-99.) # residual N fraction in roots [kg N kg-1 dry biomass]
        
        RPFLV = Float(-99.) # residual P fraction in leaves [kg P kg-1 dry biomass]
        RPFST = Float(-99.) # residual P fraction in stems [kg P kg-1 dry biomass]
        RPFRT = Float(-99.) # residual P fraction in roots [kg P kg-1 dry biomass]
        
        RKFLV = Float(-99.) # residual K fraction in leaves [kg P kg-1 dry biomass]
        RKFST = Float(-99.) # residual K fraction in stems [kg P kg-1 dry biomass]
        RKFRT = Float(-99.) # residual K fraction in roots [kg P kg-1 dry biomass]    
        
        FNTRT = Float(-99.) # NPK translocation from roots as a fraction of resp. total NPK amounts translocated from leaves and stems
        
        
    class StateVariables(StatesTemplate):
        ATNLV = Float(-99.) # translocatable N amount in leaves [kg N ha-1]
        ATNST = Float(-99.) # translocatable N amount in stems [kg N ha-1]
        ATNRT = Float(-99.) # translocatable N amount in roots [kg N ha-1]
        
        ATPLV = Float(-99.) # translocatable P amount in leaves [kg N ha-1]
        ATPST = Float(-99.) # translocatable P amount in stems [kg N ha-1]
        ATPRT = Float(-99.) # translocatable P amount in roots [kg N ha-1]
        
        ATKLV = Float(-99.) # translocatable K amount in leaves [kg N ha-1
        ATKST = Float(-99.) # translocatable K amount in stems [kg N ha-1]
        ATKRT = Float(-99.) # translocatable K amount in roots [kg N ha-1]
   
    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk translocation __call__()
        """

        self.params = self.Parameters(cropdata)
        self.states  = self.StateVariables(kiosk,publish=["ATNLV","ATNST","ATNRT",
                                                        "ATPLV","ATPST","ATPRT",
                                                        "ATKLV","ATKST","ATKRT"],
                                            ATNLV=0.,ATNST=0.,ATNRT=0.,
                                            ATPLV=0.,ATPST=0.,ATPRT=0.,
                                            ATKLV=0.,ATKST=0.,ATKRT=0.)
        self.kiosk  = kiosk
        
    @prepare_states
    def integrate(self, day):
        params = self.params
        states = self.states 
        
        WLV  = self.kiosk["WLV"]
        WST  = self.kiosk["WST"]
        WRT  = self.kiosk["WRT"]
        
        ANLV  = self.kiosk["ANLV"]
        ANST  = self.kiosk["ANST"]
        ANRT  = self.kiosk["ANRT"]
        
        APLV  = self.kiosk["APLV"]
        APST  = self.kiosk["APST"]
        APRT  = self.kiosk["APRT"]
        
        AKLV  = self.kiosk["AKLV"]
        AKST  = self.kiosk["AKST"]
        AKRT  = self.kiosk["AKRT"]
               
        
#       translocatable N amount in the organs [kg N ha-1]
        states.ATNLV = max (0. , ANLV - WLV * params.RNFLV)
        states.ATNST = max (0. , ANST - WST * params.RNFST)
        states.ATNRT = max((states.ATNLV + states.ATNST) * params.FNTRT, ANRT - WRT * params.RNFRT)
               

#       translocatable P amount in the organs [kg P ha-1]
        states.ATPLV = max (0. , APLV - WLV * params.RPFLV)
        states.ATPST = max (0. , APST - WST * params.RPFST)
        states.ATPRT = max((states.ATPLV + states.ATPST) * params.FNTRT, APRT - WRT * params.RPFRT)
        
      
#       translocatable K amount in the organs [kg K ha-1]
        states.ATKLV = max (0. , AKLV - WLV * params.RKFLV)
        states.ATKST = max (0. , AKST - WST * params.RKFST)
        states.ATKRT = max((states.ATKLV + states.ATKST) * params.FNTRT, AKRT - WRT * params.RKFRT)

        
      
     
