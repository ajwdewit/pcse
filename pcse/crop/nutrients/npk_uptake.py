
from ...traitlets import Float, Instance
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class npk_uptake_rate(SimulationObject):
    
    class Parameters(ParamTemplate):
        NFIXF  = Float(-99.) # fraction of crop nitrogen uptake by biological fixation
        DVSNT  = Float(-99.) # development stage above which NPK translocation to storage organs does occur 
        DVSNLT = Float(-99.) # development stage above which no crop N-P-K uptake does occur
    
        
    class RateVariables(RatesTemplate):
        RNULV  = Float(-99.) # N uptake rate [kg ha-1 d -1]
        RNUST  = Float(-99.)
        RNURT  = Float(-99.)
        RNUSO  = Float(-99.)
        
        RPULV  = Float(-99.) # P uptake rate [kg ha-1 d -1]
        RPUST  = Float(-99.)
        RPURT  = Float(-99.)
        RPUSO  = Float(-99.)

        RKULV  = Float(-99.) # N uptake rate [kg ha-1 d -1]
        RKUST  = Float(-99.)
        RKURT  = Float(-99.)
        RKUSO  = Float(-99.)
        
        NUPTR  = Float(-99.) # Total N uptake rate [kg ha-1 d -1]
        PUPTR  = Float(-99.)
        KUPTR  = Float(-99.)
        NFIXTR = Float(-99.)
            
    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk uptake __call__()
        """
        self.rates   = self.RateVariables(kiosk,publish=["RNULV","RNUST","RNURT","RNUSO", 
                                                         "RPULV","RPUST","RPURT","RPUSO", 
                                                         "RKULV","RKUST","RKURT","RKUSO",
                                                         "NUPTR","PUPTR","KUPTR","NFIXTR"])
        self.params = self.Parameters(cropdata)
        self.kiosk  = kiosk

    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        params = self.params
        
        NMINT  = self.kiosk["NMINT"] # total mineral N from soil and fertiliser  [kg ha-1]
        PMINT  = self.kiosk["PMINT"] # total mineral P from soil and fertiliser  [kg ha-1]
        KMINT  = self.kiosk["KMINT"] # total mineral K from soil and fertiliser  [kg ha-1]
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NDEML  = self.kiosk["NDEML"]  # N demand leaves [kg ha-1]
        NDEMS  = self.kiosk["NDEMS"]  # N demand stems [kg ha-1]
        NDEMR  = self.kiosk["NDEMR"]  # N demand roots [kg ha-1]
        NDEMSO = self.kiosk["NDEMSO"] # N demand storage organs [kg ha-1]

        PDEML  = self.kiosk["PDEML"]  # P demand leaves [kg ha-1]
        PDEMS  = self.kiosk["PDEMS"]  # P demand stems [kg ha-1]
        PDEMR  = self.kiosk["PDEMR"]  # P demand roots [kg ha-1]
        PDEMSO = self.kiosk["PDEMSO"] # P demand storage organs [kg ha-1]
 
        KDEML  = self.kiosk["KDEML"]  # K demand leaves [kg ha-1]
        KDEMS  = self.kiosk["KDEMS"]  # K demand stems [kg ha-1]
        KDEMR  = self.kiosk["KDEMR"]  # K demand roots [kg ha-1]
        KDEMSO = self.kiosk["KDEMSO"] # K demand storage organs [kg ha-1]
        
        NUPSO = self.kiosk["NUPSO"]   # N supply to storage organs [kg ha-1]
        PUPSO = self.kiosk["PUPSO"]   # P supply to storage organs [kg ha-1]
        KUPSO = self.kiosk["KUPSO"]   # K supply to storage organs [kg ha-1]
        
#       total NPK demand of leaves, stems and roots
        NDEMTO = NDEML + NDEMS + NDEMR
        PDEMTO = PDEML + PDEMS + PDEMR
        KDEMTO = KDEML + KDEMS + KDEMR  
        
#       NPK uptake rate in storage organs (kg N ha-1 d-1)
#       is the mimimum of supply and demand
        rates.RNUSO =  min(NDEMSO, NUPSO)
        rates.RPUSO =  min(PDEMSO, PUPSO)
        rates.RKUSO =  min(KDEMSO, KUPSO)
        
        TRANRF = TRA/TRAMX
        
#       No nutrients are absorbed after developmentstage DVSNLT or
#       when watershortage occurs i.e. TRANRF <= 0.01
        if DVS < params.DVSNLT and TRANRF > 0.01 :
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.

        
        # NPK uptake rate from soil
        rates.NUPTR = (max (0., min ((1.-params.NFIXF)*NDEMTO, NMINT))* NutrientLIMIT)
        rates.PUPTR = (max (0., min (PDEMTO, PMINT))* NutrientLIMIT)
        rates.KUPTR = (max (0., min (KDEMTO, KMINT))* NutrientLIMIT)
       
        
        # biological nitrogen fixation
        rates.NFIXTR= (max (0., params.NFIXF*NDEMTO)* NutrientLIMIT)
        
        # NPK uptake rate
        # if no demand then uptake rate = 0.
        if NDEMTO == 0.:
            rates.RNULV = rates.RNUST = rates.RNURT =0.
        else:    
            rates.RNULV = (NDEML / NDEMTO) * (rates.NUPTR + rates.NFIXTR)
            rates.RNUST = (NDEMS / NDEMTO) * (rates.NUPTR + rates.NFIXTR)
            rates.RNURT = (NDEMR / NDEMTO) * (rates.NUPTR + rates.NFIXTR)
            
            
        if PDEMTO == 0.:
            rates.RPULV = rates.RPUST = rates.RPURT =0.
        else:    
            rates.RPULV = (PDEML / PDEMTO) * rates.PUPTR
            rates.RPUST = (PDEMS / PDEMTO) * rates.PUPTR
            rates.RPURT = (PDEMR / PDEMTO) * rates.PUPTR
            

        if KDEMTO == 0.:
            rates.RKULV = rates.RKUST = rates.RKURT =0.
        else:    
            rates.RKULV = (KDEML / KDEMTO) * rates.KUPTR
            rates.RKUST = (KDEMS / KDEMTO) * rates.KUPTR
            rates.RKURT = (KDEMR / KDEMTO) * rates.KUPTR
            
#        print "RKUST: ",rates.RKUST, "KDEMS: ",KDEMS, "KDEMTO: ",KDEMTO, "KUPTR: ",rates.KUPTR, KMINT
        
        
        
        
        
        
        
        
        
        