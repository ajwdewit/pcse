#!/usr/bin/env python


from ...util import doy
from ...traitlets import Float, Instance, AfgenTrait
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from ... import signals

class NPK_Soil_Dynamics(SimulationObject):
    
    NSOILI = Float(-99.) # initial soil N amount
    PSOILI = Float(-99.) # initial soil P amount
    KSOILI = Float(-99.) # initial soil K amount
    
    class Parameters(ParamTemplate):      
        # FERNTAB = AfgenTrait() # N fertilizer application table
        # FERPTAB = AfgenTrait() # P fertilizer application table
        # FERKTAB = AfgenTrait() # K fertilizer application table
        #
        # NRFTAB = AfgenTrait() # kg N uptake per kg fertilser-N applied
        # PRFTAB = AfgenTrait() # kg P uptake per kg fertilser-P applied
        # KRFTAB = AfgenTrait() # kg K uptake per kg fertilser-K applied
        
        NSOILBASE = Float(-99.)  # total mineral soil N available at start of growth period [kg N/ha]
        NSOILBASE_FR = Float(-99.)  # fraction of soil mineral N coming available per day [day-1]

        PSOILBASE = Float(-99.)  # total mineral soil P available at start of growth period [kg N/ha]
        PSOILBASE_FR = Float(-99.)  # fraction of soil mineral P coming available per day [day-1]
        
        KSOILBASE = Float(-99.)  # total mineral soil K available at start of growth period [kg N/ha]
        KSOILBASE_FR = Float(-99.)  # fraction of soil mineral K coming available per day [day-1]
        
        DVSNPK_STOP = Float(-99.)  # Development stage after which no nutrients are absorbed

        # Background rates of N/P/K supply [kg/ha/day]
        BACKGROUND_N_SUPPLY = Float()
        BACKGROUND_P_SUPPLY = Float()
        BACKGROUND_K_SUPPLY = Float()

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

        # Rate of fertilizer supply for N/P/K [kg/ha/day]
        FERT_N_SUPPLY = Float()
        FERT_P_SUPPLY = Float()
        FERT_K_SUPPLY = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk
        
        # INITIAL STATES
        p = self.params
       
        self.NSOILI = p.NSOILBASE
        self.PSOILI = p.PSOILBASE
        self.KSOILI = p.KSOILBASE
        
        self.states = self.StateVariables(kiosk,
            publish=["NSOIL", "PSOIL", "KSOIL", "NAVAIL", "PAVAIL", "KAVAIL"],
            NSOIL=p.NSOILBASE, PSOIL=p.PSOILBASE, KSOIL=p.KSOILBASE,
            NAVAIL=0.091, PAVAIL=0.091, KAVAIL=0.091)

        self._connect_signal(self._on_APPLY_NPK, signals.apply_npk)
        
    @prepare_rates
    def calc_rates(self, day):
        r = self.rates
        s = self.states
        p = self.params
        
        TRA   = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        DVS   = self.kiosk["DVS"]
        
        NUPTR = self.kiosk["RNUPTAKE"]
        PUPTR = self.kiosk["RPUPTAKE"]
        KUPTR = self.kiosk["RKUPTAKE"]
        
        TRANRF = TRA/TRAMX
        
        if DVS < p.DVSNPK_STOP and TRANRF > 0.01 :
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.
                    
        r.RNSOIL = -max(0., min(p.NSOILBASE_FR * self.NSOILI * NutrientLIMIT, s.NSOIL))
        r.RPSOIL = -max(0., min(p.PSOILBASE_FR * self.PSOILI * NutrientLIMIT, s.PSOIL))
        r.RKSOIL = -max(0., min(p.KSOILBASE_FR * self.KSOILI * NutrientLIMIT, s.KSOIL))
               
        r.RNAVAIL = r.FERT_N_SUPPLY + p.BACKGROUND_N_SUPPLY - NUPTR - r.RNSOIL
        r.RPAVAIL = r.FERT_P_SUPPLY + p.BACKGROUND_P_SUPPLY - PUPTR - r.RPSOIL
        r.RKAVAIL = r.FERT_K_SUPPLY + p.BACKGROUND_K_SUPPLY - KUPTR - r.RKSOIL
        
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

    def _on_APPLY_NPK(self, N_amount=None, P_amount=None, K_amount=None, N_recovery=None,
                      P_recovery=None, K_recovery=None):

        r = self.rates
        r.unlock()
        r.FERT_K_SUPPLY = N_amount * N_recovery
        r.FERT_P_SUPPLY = P_amount * P_recovery
        r.FERT_K_SUPPLY = K_amount * K_recovery
        r.lock()