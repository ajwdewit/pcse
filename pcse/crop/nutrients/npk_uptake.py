#TODO: integrate npk_demand.py and npk_uptake.py because these are complementary
from ...traitlets import Float, Instance
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class npk_uptake_rate(SimulationObject):
    
    class Parameters(ParamTemplate):
        NFIX_FR = Float(-99.)  # fraction of crop nitrogen uptake by biological fixation
        DVSNPK_TRANSL = Float(-99.)  # development stage above which NPK translocation to storage organs does occur
        DVSNPK_STOP = Float(-99.)  # development stage above which no crop N-P-K uptake does occur

    class RateVariables(RatesTemplate):
        RNULV  = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RNUST  = Float(-99.)
        RNURT  = Float(-99.)
        RNUSO  = Float(-99.)
        
        RPULV  = Float(-99.)  # P uptake rate [kg ha-1 d -1]
        RPUST  = Float(-99.)
        RPURT  = Float(-99.)
        RPUSO  = Float(-99.)

        RKULV  = Float(-99.)  # N uptake rate [kg ha-1 d -1]
        RKUST  = Float(-99.)
        RKURT  = Float(-99.)
        RKUSO  = Float(-99.)
        
        RNUPTAKE  = Float(-99.)  # Total N uptake rate [kg ha-1 d -1]
        RPUPTAKE  = Float(-99.)
        RKUPTAKE  = Float(-99.)
        RNFIX = Float(-99.)
            
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
                                                         "RNUPTAKE","RPUPTAKE","RKUPTAKE",
                                                         "RNFIX"])
        self.params = self.Parameters(cropdata)
        self.kiosk  = kiosk

    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        params = self.params
        
        NAVAIL = self.kiosk["NAVAIL"]  # total mineral N from soil and fertiliser  [kg ha-1]
        PAVAIL = self.kiosk["PAVAIL"]  # total mineral P from soil and fertiliser  [kg ha-1]
        KAVAIL = self.kiosk["KAVAIL"]  # total mineral K from soil and fertiliser  [kg ha-1]
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NDEMLV = self.kiosk["NDEMLV"]  # N demand leaves [kg ha-1]
        NDEMST = self.kiosk["NDEMST"]  # N demand stems [kg ha-1]
        NDEMRT = self.kiosk["NDEMRT"]  # N demand roots [kg ha-1]
        NDEMSO = self.kiosk["NDEMSO"]  # N demand storage organs [kg ha-1]

        PDEMLV = self.kiosk["PDEMLV"]  # P demand leaves [kg ha-1]
        PDEMST = self.kiosk["PDEMST"]  # P demand stems [kg ha-1]
        PDEMRT = self.kiosk["PDEMRT"]  # P demand roots [kg ha-1]
        PDEMSO = self.kiosk["PDEMSO"]  # P demand storage organs [kg ha-1]
 
        KDEMLV = self.kiosk["KDEMLV"]  # K demand leaves [kg ha-1]
        KDEMST = self.kiosk["KDEMST"]  # K demand stems [kg ha-1]
        KDEMRT = self.kiosk["KDEMRT"]  # K demand roots [kg ha-1]
        KDEMSO = self.kiosk["KDEMSO"]  # K demand storage organs [kg ha-1]
        
        NTRANSLOCATABLE = self.kiosk["NTRANSLOCATABLE"]  # N supply to storage organs [kg ha-1]
        PTRANSLOCATABLE = self.kiosk["PTRANSLOCATABLE"]  # P supply to storage organs [kg ha-1]
        KTRANSLOCATABLE = self.kiosk["KTRANSLOCATABLE"]  # K supply to storage organs [kg ha-1]
        
#       total NPK demand of leaves, stems and roots
        NDEMTO = NDEMLV + NDEMST + NDEMRT
        PDEMTO = PDEMLV + PDEMST + PDEMRT
        KDEMTO = KDEMLV + KDEMST + KDEMRT
        
#       NPK uptake rate in storage organs (kg N ha-1 d-1)
#       is the mimimum of supply and demand
        rates.RNUSO = min(NDEMSO, NTRANSLOCATABLE)
        rates.RPUSO = min(PDEMSO, PTRANSLOCATABLE)
        rates.RKUSO = min(KDEMSO, KTRANSLOCATABLE)
        
        TRANRF = TRA/TRAMX
        
#       No nutrients are absorbed after developmentstage DVSNLT or
#       when watershortage occurs i.e. TRANRF <= 0.01
        if DVS < params.DVSNPK_STOP and TRANRF > 0.01:
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.

        
        # NPK uptake rate from soil
        rates.RNUPTAKE = (max (0., min((1. - params.NFIX_FR)*NDEMTO, NAVAIL)) * NutrientLIMIT)
        rates.RPUPTAKE = (max (0., min(PDEMTO, PAVAIL)) * NutrientLIMIT)
        rates.RKUPTAKE = (max (0., min(KDEMTO, KAVAIL)) * NutrientLIMIT)
       
        
        # biological nitrogen fixation
        rates.RNFIX = (max(0., params.NFIX_FR * NDEMTO) * NutrientLIMIT)
        
        # NPK uptake rate
        # if no demand then uptake rate = 0.
        if NDEMTO == 0.:
            rates.RNULV = rates.RNUST = rates.RNURT = 0.
        else:    
            rates.RNULV = (NDEMLV / NDEMTO) * (rates.RNUPTAKE + rates.RNFIX)
            rates.RNUST = (NDEMST / NDEMTO) * (rates.RNUPTAKE + rates.RNFIX)
            rates.RNURT = (NDEMRT / NDEMTO) * (rates.RNUPTAKE + rates.RNFIX)
            
            
        if PDEMTO == 0.:
            rates.RPULV = rates.RPUST = rates.RPURT = 0.
        else:    
            rates.RPULV = (PDEMLV / PDEMTO) * rates.RPUPTAKE
            rates.RPUST = (PDEMST / PDEMTO) * rates.RPUPTAKE
            rates.RPURT = (PDEMRT / PDEMTO) * rates.RPUPTAKE
            

        if KDEMTO == 0.:
            rates.RKULV = rates.RKUST = rates.RKURT = 0.
        else:    
            rates.RKULV = (KDEMLV / KDEMTO) * rates.RKUPTAKE
            rates.RKUST = (KDEMST / KDEMTO) * rates.RKUPTAKE
            rates.RKURT = (KDEMRT / KDEMTO) * rates.RKUPTAKE
            
#        print "RKUST: ",rates.RKUST, "KDEMS: ",KDEMS, "KDEMTO: ",KDEMTO, "KUPTR: ",rates.KUPTR, KMINT
