#!/usr/bin/env python


from ...util import doy
from ...traitlets import Float, Instance, AfgenTrait
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from copy import deepcopy

class NPK_Soil_Dynamics(SimulationObject):
    
    NSOILI = Float(-99.) # initial soil N amount
    PSOILI = Float(-99.) # initial soil P amount
    KSOILI = Float(-99.) # initial soil K amount
    
    class Parameters(ParamTemplate):      
        FERNTAB = AfgenTrait() # N fertilizer application table
        FERPTAB = AfgenTrait() # P fertilizer application table
        FERKTAB = AfgenTrait() # K fertilizer application table
        
        NRFTAB = AfgenTrait() # kg N uptake per kg fertilser-N applied
        PRFTAB = AfgenTrait() # kg P uptake per kg fertilser-P applied
        KRFTAB = AfgenTrait() # kg K uptake per kg fertilser-K applied
        
        NSOILBASE = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        NSOILBASE_FR = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]

        PSOILBASE = Float(-99.)  # total mineral soil P available at start of growth period [kg N/ha]
        PSOILBASE_FR = Float(-99.)  # fraction of soil mineral P coming available per day [day-1]
        
        KSOILBASE = Float(-99.)  # total mineral soil K available at start of growth period [kg N/ha]
        KSOILBASE_FR = Float(-99.)  # fraction of soil mineral K coming available per day [day-1]
        
        DVSNPK_STOP = Float(-99.)  # Developmentstage after which no nutrients are absorbed
        
    class StateVariables(StatesTemplate):
        NSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        PSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        KSOIL = Float(-99.)  # mineral N available from soil for crop    kg N ha-1
        
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1
        PAVAIL = Float(-99.)  # total mineral P from soil and fertiliser  kg N ha-1
        KAVAIL = Float(-99.)  # total mineral K from soil and fertiliser  kg N ha-1
      
    class RateVariables(RatesTemplate):
        RNSOIL = Float(-99.)
        RPSOIL = Float(-99.)
        RKSOIL = Float(-99.)
        
        RNAVAIL = Float(-99.)
        RPAVAIL = Float(-99.)
        RKAVAIL = Float(-99.)
        
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
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
      # INITIAL STATES
        params = self.params
       
        iday  = doy(day)
        NAVAIL = params.FERNTAB(iday) * params.NRFTAB(iday)
        PAVAIL = params.FERPTAB(iday) * params.PRFTAB(iday)
        KAVAIL = params.FERKTAB(iday) * params.KRFTAB(iday)
        
        self.NSOILI = params.NSOILBASE
        self.PSOILI = params.PSOILBASE
        self.KSOILI = params.KSOILBASE
        
        self.states = self.StateVariables(kiosk,
            publish=["NSOIL", "PSOIL", "KSOIL", "NAVAIL", "PAVAIL", "KAVAIL"],
            NSOIL=params.NSOILBASE, PSOIL=params.PSOILBASE, KSOIL=params.KSOILBASE,
            NAVAIL=NAVAIL, PAVAIL=PAVAIL, KAVAIL=KAVAIL)
        
        
    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        states = self.states
        params = self.params
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NUPTR = self.kiosk["RNUPTAKE"]
        PUPTR = self.kiosk["RPUPTAKE"]
        KUPTR = self.kiosk["RKUPTAKE"]
        
        TRANRF = TRA/TRAMX
        
        if DVS < params.DVSNPK_STOP and TRANRF > 0.01 :
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.
                    
        rates.RNSOIL = -max(0., min(params.NSOILBASE_FR * self.NSOILI * NutrientLIMIT, states.NSOIL))
        rates.RPSOIL = -max(0., min(params.PSOILBASE_FR * self.PSOILI * NutrientLIMIT, states.PSOIL))
        rates.RKSOIL = -max(0., min(params.KSOILBASE_FR * self.KSOILI * NutrientLIMIT, states.KSOIL))
        
        iday  = doy(day)
        FERTNS = params.FERNTAB(iday) * params.NRFTAB(iday)
        FERTPS = params.FERPTAB(iday) * params.PRFTAB(iday)
        FERTKS = params.FERKTAB(iday) * params.KRFTAB(iday)
#        print iday,FERTNS, FERTPS, FERTKS
      
               
        rates.RNAVAIL = FERTNS - NUPTR - rates.RNSOIL
        rates.RPAVAIL = FERTPS - PUPTR - rates.RPSOIL
        rates.RKAVAIL = FERTKS - KUPTR - rates.RKSOIL
        
    @prepare_states
    def integrate(self, day):
        rates = self.rates
        states = self.states

        # mineral NPK  amount in the soil
        states.NSOIL += rates.RNSOIL
        states.PSOIL += rates.RPSOIL
        states.KSOIL += rates.RKSOIL
        
        # total (soil + fertilizer) NPK amount in soil
        states.NAVAIL += rates.RNAVAIL
        states.PAVAIL += rates.RPAVAIL
        states.KAVAIL += rates.RKAVAIL
        