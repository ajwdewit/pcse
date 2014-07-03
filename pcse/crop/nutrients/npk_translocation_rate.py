#!/usr/bin/env python
"""
Class to estimate the nutrient allocation rates from the plant organs to the
storage organs
"""


from ...decorators import prepare_rates, prepare_states
from ...traitlets import Float, Instance
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class npk_translocation(SimulationObject):
    
    class RateVariables(RatesTemplate):
        RNTLV = Float(-99.) # N translocation rate from leaves [kg ha-1 d-1]
        RNTST = Float(-99.) # N translocation rate from stems [kg ha-1 d-1]
        RNTRT = Float(-99.) # N translocation rate from roots [kg ha-1 d-1]
        
        RPTLV = Float(-99.) # P translocation rate from leaves [kg ha-1 d-1]
        RPTST = Float(-99.) # P translocation rate from stems [kg ha-1 d-1]
        RPTRT = Float(-99.) # P translocation rate from roots [kg ha-1 d-1]
        
        RKTLV = Float(-99.) # K translocation rate from leaves [kg ha-1 d-1]
        RKTST = Float(-99.) # K translocation rate from stems [kg ha-1 d-1]
        RKTRT = Float(-99.) # K translocation rate from roots [kg ha-1 d-1]
   
    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk translocation __call__()
        """
        
        self.rates   = self.RateVariables(kiosk,publish=["RNTLV","RNTST","RNTRT",
                                                         "RPTLV","RPTST","RPTRT",
                                                         "RKTLV","RKTST","RKTRT"])
        self.kiosk  = kiosk
        
    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        params = self.params
        
        WLV  = self.kiosk["WLV"]
        WST  = self.kiosk["WST"]
        WRT  = self.kiosk["WRT"]
        
        ATNLV  = self.kiosk["ATNLV"] # translocatble N amount 
        ATNST  = self.kiosk["ATNST"]
        ATNRT  = self.kiosk["ATNRT"]
        
        ATPLV  = self.kiosk["ATPLV"] # translocatble P amount
        ATPST  = self.kiosk["ATPST"]
        ATPRT  = self.kiosk["ATPRT"]
        
        ATKLV  = self.kiosk["ATKLV"] # translocatble K amount
        ATKST  = self.kiosk["ATKST"]
        ATKRT  = self.kiosk["ATKRT"]
        
        RNUSO = self.kiosk["RNUSO"] # N uptake storage organs
        RPUSO = self.kiosk["RPUSO"] # P uptake storage organs
        RKUSO = self.kiosk["RKUSO"] # K uptake storage organs
               
        
#       max amount translocatable NPK [kg ha-1 d-1]
        ATN   = ATNLV + ATNST + ATNRT
        ATP   = ATPLV + ATPST + ATPRT
        ATK   = ATKLV + ATKST + ATKRT

  
#       partionioning of the uptake      
        # if max amount of translocatable N = 0 then
        # N translocation rate is 0
        if ATN > 0.:
            rates.RNTLV = RNUSO * ATNLV / ATN
            rates.RNTST = RNUSO * ATNST / ATN
            rates.RNTRT = RNUSO * ATNRT / ATN
        else:
            rates.RNTLV = rates.RNTST = rates.RNTRT = 0.

        # if max amount of translocatable P = 0 then
        # P translocation rate is 0            
        if ATP > 0:
            rates.RPTLV = RPUSO * ATPLV / ATP
            rates.RPTST = RPUSO * ATPST / ATP
            rates.RPTRT = RPUSO * ATPRT / ATP
        else:
            rates.RPTLV =  rates.RPTST = rates.RPTRT = 0.
           
        # if max amount of translocatable K = 0 then
        # K translocation rate is 0           
        if ATK > 0:
            rates.RKTLV = RKUSO * ATKLV / ATK
            rates.RKTST = RKUSO * ATKST / ATK
            rates.RKTRT = RKUSO * ATKRT / ATK
        else:
            rates.RKTLV = rates.RKTST = rates.RKTRT = 0.
            
        
