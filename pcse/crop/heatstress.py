# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2016
# Approach based on LINTUL4-vsht made by Joost Wolf

"""Implementation of a models for heatstress around flowering in WOFOST

Classes defined here:
- HeatStressFlowering:
"""
import datetime

from ..traitlets import Float, Int, Instance, Enum, Bool, AfgenTrait
from ..decorators import prepare_rates, prepare_states

from ..util import limit, daylength
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject, VariableKiosk
from .. import signals
from .. import exceptions as exc

#-------------------------------------------------------------------------------
class HeatStress_Around_Flowering(SimulationObject):
    """Implements the algorithms for heatstress around flowring in WOFOST.
    



    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    DTEMP    Average daytime temperature                    SCr        |C| day
    TBASEM   Base temperature for emergence                 SCr        |C|
    TEFFMX   Maximum effective temperature for emergence    SCr        |C|
    TSUM1    Temperature sum from emergence to anthesis     SCr        |C| day
    TSUM2    Temperature sum from anthesis to maturity      SCr        |C| day
    IDSL     Switch for phenological development options    SCr        -
             temperature only (IDSL=0), including           SCr
             daylength (IDSL=1) and including
             vernalization (IDSL>=2)
    DLO      Optimal daylength for phenological             SCr        hr
             development
    DLC      Critical daylength for phenological            SCr        hr
             development
    DVSI     Initial development stage at emergence.        SCr        -
             Usually this is zero, but it can be higher
             for crops that are transplanted (e.g. paddy
             rice)
    DVSEND   Final development stage                        SCr        -
    DTSMTB   Daily increase in temperature sum as a         TCr        |C|
             function of daily mean temperature.
    =======  ============================================= =======  ============

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    DVS      Development stage                                  Y    - 
    TSUM     Temperature sum                                    N    |C| day
    TSUME    Temperature sum for emergence                      N    |C| day
    DOS      Day of sowing                                      N    - 
    DOE      Day of emergence                                   N    - 
    DOA      Day of Anthesis                                    N    - 
    DOM      Day of maturity                                    N    - 
    DOH      Day of harvest                                     N    -
    STAGE    Current phenological stage, can take the           N    -
             folowing values:
             `emerging|vegetative|reproductive|mature`
    =======  ================================================= ==== ============

    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    DTSUME   Increase in temperature sum for emergence          N    |C|
    DTSUM    Increase in temperature sum for anthesis or        N    |C|
             maturity
    DVR      Development rate                                   Y    |day-1|
    =======  ================================================= ==== ============
    
    **External dependencies:**

    None    

    **Signals sent or handled**
    
    `DVS_Phenology` sends the `crop_finish` signal when maturity is
    reached and the `end_type` is 'maturity' or 'earliest'.
    
    """

    class Parameters(ParamTemplate):
        DVSHEB = Float(-99.)  # beginning of the period with possible heat stress effects
        DVSHEF = Float(-99.)  # finalization of the period with possible heat stress effects
        RDGRTB = AfgenTrait() # Temperature response function for grain formation as a function of day-time temp.

    #-------------------------------------------------------------------------------
    class RateVariables(RatesTemplate):
        DAYTEMP = Float(-99.)  # increase in day-time temperature sum


    #-------------------------------------------------------------------------------
    class StateVariables(StatesTemplate):
        CMDTEMP   = Float(-99.)  # Sum of the day-time temperature
        RDGRHT    = Float(-99.)   # reduction factor of heat stress round flowering

        # States which register phenological events
        DSB = Instance(datetime.date) # Day of sensitivity begin
        DSE = Instance(datetime.date) # Day of sensitivity end

        HSTAGE = Enum([None, "not sensitive", "heat sensitive", "end heat sensitivty"])

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk)
        self.kiosk  = kiosk

        #self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

        DSB, DSE, HSTAGE, RDGRHT = self._get_initial_stage(day)
        self.states = self.StateVariables(kiosk, publish="RDGRHT",RDGRHT=RDGRHT,
                                          CMDTEMP=0., DSB=DSB, DSE=DSE, HSTAGE=HSTAGE)

    #---------------------------------------------------------------------------
    def _get_initial_stage(self, day):
        """"""
        p = self.params

        HSTAGE = "not sensitive"
        DSB = None
        DSE = None
        RDGRHT = 1.
            
        return (DSB, DSE, HSTAGE, RDGRHT)

    @prepare_rates
    def calc_rates(self, day, drv):
        r     = self.rates
        kiosk = self.kiosk

        r.DAYTEMP = drv.DTEMP
        #self.rates = self.RateVariables(kiosk,DAYTEMP=DAYTEMP)

    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day):
        """Updates the state variable and checks for phenologic stages
        """

        p = self.params
        r = self.rates
        s = self.states

        DVS = self.kiosk["DVS"]

        if s.HSTAGE == "not sensitive":
            s.CMDTEMP = 0.
            if DVS >= p.DVSHEB:
                self._next_stage(day)
                RDGHRT = 1.

        elif s.HSTAGE == "heat sensitive":
            s.CMDTEMP += r.DAYTEMP
            if DVS >= p.DVSHEF:
                self._next_stage(day)
                s.RDGRHT = self._get_heatstress_factor(s.CMDTEMP)
                
        elif s.HSTAGE == "end heat sensitivty":
            s.RDGHRT = s.RDGRHT

        else: # Problem no stage defined
            msg = "No HSTAGE defined in the heat stress around flowering module."
            raise exc.PCSEError(msg)
            
        msg = "Finished state integration for %s"
        self.logger.debug(msg % day)

    #---------------------------------------------------------------------------
    def _next_stage(self, day):
        """Moves states.STAGE to the next heat sensitivity stage"""
        s = self.states

        current_HSTAGE = s.HSTAGE
        if s.HSTAGE == "not sensitive":
            s.HSTAGE = "heat sensitive"
            s.DSB = day
            
        elif s.HSTAGE == "heat sensitive":
            s.HSTAGE = "end heat sensitivty"
            s.DSE = day

        else: # Problem no stage defined
            msg = "No HSTAGE defined in the heat stress around flowering module."
            raise exc.PCSEError(msg)
        
        msg = "Changed heat sensivity stage '%s' to '%s' on %s"
        self.logger.info(msg % (current_HSTAGE, s.HSTAGE, day))

    #---------------------------------------------------------------------------
    def _get_heatstress_factor(self, day):
        p = self.params
        s = self.states

        CDAYS  = s.DSE - s.DSB
        MDTEMP = s.CMDTEMP / max(1., CDAYS)

        return p.RDGRTB(MDTEMP)



