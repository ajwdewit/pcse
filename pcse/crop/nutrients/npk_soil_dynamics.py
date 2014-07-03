#!/usr/bin/env python


from ...util import doy
from ...traitlets import Float, Instance, AfgenTrait
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from copy import deepcopy

class NPK_Soil_Dynamics(SimulationObject):
    
    NMINI = Float(-99.) # initial soil N amount
    PMINI = Float(-99.) # initial soil P amount
    KMINI = Float(-99.) # initial soil K amount
    
    class Parameters(ParamTemplate):      
        FERNTAB = AfgenTrait() # N fertilizer application table
        FERPTAB = AfgenTrait() # P fertilizer application table
        FERKTAB = AfgenTrait() # K fertilizer application table
        
        NRFTAB  = AfgenTrait() # kg N uptake per kg fertilser-N applied
        PRFTAB  = AfgenTrait() # kg P uptake per kg fertilser-P applied
        KRFTAB  = AfgenTrait() # kg K uptake per kg fertilser-K applied
        
        NMINS   = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        RTNMINS = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]

        PMINS   = Float(-99.)  # total mineral soil P available at start of growth period [kg N/ha]
        RTPMINS = Float(-99.)  # fraction of soil mineral P coming available per day [day-1]
        
        KMINS   = Float(-99.)  # total mineral soil K available at start of growth period [kg N/ha]
        RTKMINS = Float(-99.)  # fraction of soil mineral K coming available per day [day-1]
        
        DVSNLT  = Float(-99.)  # Developmentstage after which no nutrients are absorbed
        
    class StateVariables(StatesTemplate):
        NMIN = Float(-99.) # mineral N available from soil for crop    kg N ha-1
        PMIN = Float(-99.) # mineral N available from soil for crop    kg N ha-1
        KMIN = Float(-99.) # mineral N available from soil for crop    kg N ha-1
        
        NMINT = Float(-99.) # total mineral N from soil and fertiliser  kg N ha-1
        PMINT = Float(-99.) # total mineral P from soil and fertiliser  kg N ha-1
        KMINT = Float(-99.) # total mineral K from soil and fertiliser  kg N ha-1
      
    class RateVariables(RatesTemplate):
        RNMINS = Float(-99.)
        RPMINS = Float(-99.)
        RKMINS = Float(-99.)
        
        RNMINT = Float(-99.)
        RPMINT = Float(-99.)
        RKMINT = Float(-99.)
        
    def initialize(self, day, kiosk, cropdata, fertilizer):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """
        # Merge dictionaries in order to pass them to the Parameters class.
        # use merge_dict iso deepcopy?
        parvalues = deepcopy(cropdata)
        parvalues.update(fertilizer)
        
        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk)
        self.kiosk  = kiosk
        
      # INITIAL STATES
        params= self.params
       
        NMIN = params.NMINS
        PMIN = params.PMINS
        KMIN = params.KMINS
        
        iday  = doy(day)
        NMINT = params.FERNTAB(iday) * params.NRFTAB(iday)
        PMINT = params.FERPTAB(iday) * params.PRFTAB(iday)
        KMINT = params.FERKTAB(iday) * params.KRFTAB(iday)
        
        self.NMINI = NMIN
        self.PMINI = PMIN
        self.KMINI = KMIN
        
        self.states = self.StateVariables(kiosk, publish=["NMIN","PMIN","KMIN",
                                                          "NMINT","PMINT","KMINT"],
            NMIN=NMIN, PMIN=PMIN, KMIN=KMIN,NMINT=NMINT, PMINT=PMINT, KMINT=KMINT)
        
        
    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        states = self.states
        params = self.params
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NMIN  = self.kiosk["NMIN"]
        PMIN  = self.kiosk["PMIN"]
        KMIN  = self.kiosk["KMIN"]
        
        NUPTR = self.kiosk["NUPTR"]
        PUPTR = self.kiosk["PUPTR"]
        KUPTR = self.kiosk["KUPTR"]
        
        TRANRF = TRA/TRAMX
        
        if DVS < params.DVSNLT and TRANRF > 0.01 :
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.
                    
        rates.RNMINS = -max(0.,min( params.RTNMINS * self.NMINI * NutrientLIMIT, NMIN))
        rates.RPMINS = -max(0.,min( params.RTPMINS * self.PMINI * NutrientLIMIT, PMIN))
        rates.RKMINS = -max(0.,min( params.RTKMINS * self.KMINI * NutrientLIMIT, KMIN))
        
        iday  = doy(day)
        FERTNS = params.FERNTAB(iday) * params.NRFTAB(iday)
        FERTPS = params.FERPTAB(iday) * params.PRFTAB(iday)
        FERTKS = params.FERKTAB(iday) * params.KRFTAB(iday)
#        print iday,FERTNS, FERTPS, FERTKS
      
               
        rates.RNMINT = FERTNS - NUPTR - rates.RNMINS
        rates.RPMINT = FERTPS - PUPTR - rates.RPMINS
        rates.RKMINT = FERTKS - KUPTR - rates.RKMINS
        
    @prepare_states
    def integrate(self, day):
        rates  = self.rates
        states = self.states

        # mineral NPK  amount in the soil
        states.NMIN += rates.RNMINS
        states.PMIN += rates.RPMINS
        states.KMIN += rates.RKMINS
        
        # total (soil + fertilizer) NPK amount in soil
        states.NMINT += rates.RNMINT
        states.PMINT += rates.RPMINT
        states.KMINT += rates.RKMINT
        