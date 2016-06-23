# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Implementation of a models for phenological development in WOFOST

Classes defined here:
- DVS_Phenology: Implements the algorithms for phenologic development
- Vernalisation: 
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
class Vernalisation(SimulationObject):
    """ Modification of phenological development due to vernalisation.
    
    The vernalization approach here is based on the work of Lenny van Bussel
    (2011), which in turn is based on Wang and Engel (1998). The basic
    principle is that winter wheat needs a certain number of days with temperatures
    within an optimum temperature range to complete its vernalisation
    requirement. Until the vernalisation requirement is fulfilled, the crop
    development is delayed.
    
    The rate of vernalization (VERNR) is defined by the temperature response
    function VERNRTB. Within the optimal temperature range 1 day is added
    to the vernalisation state (VERN). The reduction on the phenological
    development is calculated from the base and saturated vernalisation
    requirements (VERNBASE and VERNSAT). The reduction factor (VERNFAC) is
    scaled linearly between VERNBASE and VERNSAT.
    
    A critical development stage (VERNDVS) is used to stop the effect of
    vernalisation when this DVS is reached. This is done to improve model
    stability in order to avoid that Anthesis is never reached due to a
    somewhat too high VERNSAT. Nevertheless, a warning is written to the log
    file, if this happens.   
    
    * Van Bussel, 2011. From field to globe: Upscaling of crop growth modelling.
      Wageningen PhD thesis. http://edepot.wur.nl/180295
    * Wang and Engel, 1998. Simulation of phenological development of wheat
      crops. Agric. Systems 58:1 pp 1-24

    *Simulation parameters* (provide in cropdata dictionary)
    
    ======== ============================================= =======  ============
     Name     Description                                   Type     Unit
    ======== ============================================= =======  ============
    VERNSAT  Saturated vernalisation requirements           SCr        days
    VERNBASE Base vernalisation requirements                SCr        days
    VERNRTB  Rate of vernalisation as a function of daily   TCr        -
             mean temperature.
    VERNDVS  Critical development stage after which the     SCr        -
             effect of vernalisation is halted
    ======== ============================================= =======  ============

    **State variables**

    ============ ================================================= ==== ========
     Name        Description                                       Pbl   Unit
    ============ ================================================= ==== ========
    VERN         Vernalisation state                                N    days
    DOV          Day when vernalisation requirements are            N    -
                 fulfilled.
    ISVERNALISED Flag indicated that vernalisation                  Y    -
                 requirement has been reached
    ============ ================================================= ==== ========


    **Rate variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    VERNR    Rate of vernalisation                              N     -
    VERNFAC  Reduction factor on development rate due to        Y     -
             vernalisation effect.
    =======  ================================================= ==== ============

    
    **External dependencies:**

    ============ =============================== ========================== =====
     Name        Description                         Provided by             Unit
    ============ =============================== ========================== =====
    DVS          Development Stage                 Phenology                 -
                 Used only to determine if the
                 critical development stage for
                 vernalisation (VERNDVS) is
                 reached.
    ============ =============================== ========================== =====
    """
    # Helper variable to indicate that DVS > VERNDVS
    _force_vernalisation = Bool(False)

    class Parameters(ParamTemplate):
        VERNSAT = Float(-99.)     # Saturated vernalisation requirements
        VERNBASE = Float(-99.)     # Base vernalisation requirements
        VERNRTB = AfgenTrait()    # Vernalisation temperature response
        VERNDVS = Float(-99.)     # Critical DVS for vernalisation fulfillment

    class RateVariables(RatesTemplate):
        VERNR = Float(-99.)        # Rate of vernalisation
        VERNFAC = Float(-99.)      # Red. factor for phenol. devel.

    class StateVariables(StatesTemplate):
        VERN = Float(-99.)              # Vernalisation state
        DOV = Instance(datetime.date)  # Day when vernalisation
                                            # requirements are fulfilled
        ISVERNALISED =  Bool()              # True when VERNSAT is reached and
                                            # Forced when DVS > VERNDVS

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs

        """
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["VERNFAC"])
        self.kiosk = kiosk

        # Define initial states
        self.states = self.StateVariables(kiosk, VERN=0., VERNFAC=0.,
                                          DOV=None, ISVERNALISED=False,
                                          publish=["ISVERNALISED"])
    #---------------------------------------------------------------------------
    @prepare_rates
    def calc_rates(self, day, drv):
        rates = self.rates
        states = self.states
        params = self.params

        DVS = self.kiosk["DVS"]
        if not states.ISVERNALISED:
            if DVS < params.VERNDVS:
                rates.VERNR = params.VERNRTB(drv.TEMP)
                r = (states.VERN - params.VERNBASE)/(params.VERNSAT-params.VERNBASE)
                rates.VERNFAC = limit(0., 1., r)
            else:
                rates.VERNR = 0.
                rates.VERNFAC = 1.0
                self._force_vernalisation = True
        else:
            rates.VERNR = 0.
            rates.VERNFAC = 1.0
    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day, delt=1.0):
        states = self.states
        rates = self.rates
        params = self.params
        
        states.VERN += rates.VERNR
        
        if states.VERN >= params.VERNSAT:  # Vernalisation requirements reached
            states.ISVERNALISED = True
            if states.DOV is None:
                states.DOV = day
                msg = "Vernalization requirements reached at day %s."
                self.logger.info(msg % day)

        elif self._force_vernalisation:  # Critical DVS for vernalisation reached
            # Force vernalisation, but do not set DOV
            states.ISVERNALISED = True

            # Write log message to warn about forced vernalisation
            msg = ("Critical DVS for vernalization (VERNDVS) reached " +
                   "at day %s, " +
                   "but vernalization requirements not yet fulfilled. " +
                   "Forcing vernalization now (VERN=%f).")
            self.logger.warning(msg % (day, states.VERN))

        else:  # Reduction factor for phenologic development
            states.ISVERNALISED = False

#-------------------------------------------------------------------------------
class DVS_Phenology(SimulationObject):
    """Implements the algorithms for phenologic development in WOFOST.
    
    Phenologic development in WOFOST is expresses using a unitless scale which
    takes the values 0 at emergence, 1 at Anthesis (flowering) and 2 at
    maturity. This type of phenological development is mainly representative
    for cereal crops. All other crops that are simulated with WOFOST are
    forced into this scheme as well, although this may not be appropriate for
    all crops. For example, for potatoes development stage 1 represents the
    start of tuber formation rather than flowering.
    
    Phenological development is mainly governed by temperature and can be
    modified by the effects of day length and vernalization 
    during the period before Anthesis. After Anthesis, only temperature
    influences the development rate.


    **Simulation parameters**
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    TSUMEM   Temperature sum from sowing to emergence       SCr        |C| day
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
    # Placeholder for start/stop types and vernalisation module
    vernalisation = Instance(Vernalisation)

    class Parameters(ParamTemplate):
        TSUMEM = Float(-99.)  # Temp. sum for emergence
        TBASEM = Float(-99.)  # Base temp. for emergence
        TEFFMX = Float(-99.)  # Max eff temperature for emergence
        TSUM1  = Float(-99.)  # Temperature sum emergence to anthesis
        TSUM2  = Float(-99.)  # Temperature sum anthesis to maturity
        IDSL   = Float(-99.)  # Switch for photoperiod (1) and vernalisation (2)
        DLO    = Float(-99.)  # Optimal day length for phenol. development
        DLC    = Float(-99.)  # Critical day length for phenol. development
        DVSI   = Float(-99.)  # Initial development stage
        DVSEND = Float(-99.)  # Final development stage
        DTSMTB = AfgenTrait() # Temperature response function for phenol.
                              # development.
        CROP_START_TYPE = Enum(["sowing", "emergence"])
        CROP_END_TYPE = Enum(["maturity", "harvest", "earliest"])

    #-------------------------------------------------------------------------------
    class RateVariables(RatesTemplate):
        DTSUME = Float(-99.)  # increase in temperature sum for emergence
        DTSUM  = Float(-99.)  # increase in temperature sum
        DVR    = Float(-99.)  # development rate

    #-------------------------------------------------------------------------------
    class StateVariables(StatesTemplate):
        DVS   = Float(-99.)  # Development stage
        TSUM  = Float(-99.)  # Temperature sum state
        TSUME = Float(-99.)  # Temperature sum for emergence state
        # States which register phenological events
        DOS = Instance(datetime.date) # Day of sowing
        DOE = Instance(datetime.date) # Day of emergence
        DOA = Instance(datetime.date) # Day of anthesis
        DOM = Instance(datetime.date) # Day of maturity
        DOH = Instance(datetime.date) # Day of harvest
        STAGE  = Enum([None, "emerging", "vegetative", "reproductive", "mature"])

    #---------------------------------------------------------------------------
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.kiosk = kiosk

        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

        # Define initial states
        DVS = self.params.DVSI
        DOS, DOE, STAGE = self._get_initial_stage(day)
        self.states = self.StateVariables(kiosk, publish="DVS",
                                          TSUM=0., TSUME=0., DVS=DVS,
                                          DOS=DOS, DOE=DOE, DOA=None, DOM=None,
                                          DOH=None, STAGE=STAGE)

        # initialize vernalisation for IDSL=2
        if self.params.IDSL >= 2:
            self.vernalisation = Vernalisation(day, kiosk, parvalues)
    
    #---------------------------------------------------------------------------
    def _get_initial_stage(self, day):
        """"""
        p = self.params

        # Define initial stage type (emergence/sowing) and fill the
        # respective day of sowing/emergence (DOS/DOE)
        if p.CROP_START_TYPE == "emergence":
            STAGE = "vegetative"
            DOE = day
            DOS = None

            # send signal to indicate crop emergence
            self._send_signal(signals.crop_emerged)

        elif p.CROP_START_TYPE == "sowing":
            STAGE = "emerging"
            DOS = day
            DOE = None

        else:
            msg = "Unknown start type: %s" % p.CROP_START_TYPE
            raise exc.PCSEError(msg)
            
        return (DOS, DOE, STAGE)

    #---------------------------------------------------------------------------
    @prepare_rates
    def calc_rates(self, day, drv):
        """Calculates the rates for phenological development
        """
        p = self.params
        r = self.rates
        s = self.states

        # Day length sensitivity
        DVRED = 1.
        if p.IDSL >= 1:
            DAYLP = daylength(day, drv.LAT)
            DVRED = limit(0., 1., (DAYLP - p.DLC)/(p.DLO - p.DLC))

        # Vernalisation
        VERNFAC = 1.
        if p.IDSL >= 2:
            if s.STAGE == 'vegetative':
                self.vernalisation.calc_rates(day, drv)
                VERNFAC = self.kiosk["VERNFAC"]

        # Development rates
        if s.STAGE == "emerging":
            r.DTSUME = limit(0., (p.TEFFMX - p.TBASEM), (drv.TEMP - p.TBASEM))
            r.DTSUM = 0.
            r.DVR = 0.

        elif s.STAGE == 'vegetative':
            r.DTSUME = 0.
            r.DTSUM = p.DTSMTB(drv.TEMP) * VERNFAC * DVRED
            r.DVR = r.DTSUM/p.TSUM1

        elif s.STAGE in ['reproductive', 'mature']:
            r.DTSUME = 0.
            r.DTSUM = p.DTSMTB(drv.TEMP)
            r.DVR = r.DTSUM/p.TSUM2

        else: # Problem: no stage defined
            msg = "Unrecognized STAGE defined in phenology submodule: %s"
            raise exc.PCSEError(msg, self.states.STAGE)
        
        msg = "Finished rate calculation for %s"
        self.logger.debug(msg % day)
        
    #---------------------------------------------------------------------------
    @prepare_states
    def integrate(self, day, delt=1.0):
        """Updates the state variable and checks for phenologic stages
        """

        p = self.params
        r = self.rates
        s = self.states

        # Integrate vernalisation module
        if p.IDSL >= 2:
            if s.STAGE == 'vegetative':
                self.vernalisation.integrate(day, delt)
            else:
                self.vernalisation.touch()

        # Integrate phenologic states
        s.TSUME += r.DTSUME
        s.DVS += r.DVR
        s.TSUM += r.DTSUM

        # Check if a new stage is reached
        if s.STAGE == "emerging":
            if s.TSUME >= p.TSUMEM:
                self._next_stage(day)
        elif s.STAGE == 'vegetative':
            if s.DVS >= 1.0:
                self._next_stage(day)
                s.DVS = 1.0
        elif s.STAGE == 'reproductive':
            if s.DVS >= p.DVSEND:
                self._next_stage(day)
        elif s.STAGE == 'mature':
            pass
        else: # Problem no stage defined
            msg = "No STAGE defined in phenology submodule"
            raise exc.PCSEError(msg)
            
        msg = "Finished state integration for %s"
        self.logger.debug(msg % day)

    #---------------------------------------------------------------------------
    def _next_stage(self, day):
        """Moves states.STAGE to the next phenological stage"""
        s = self.states
        p = self.params

        current_STAGE = s.STAGE
        if s.STAGE == "emerging":
            s.STAGE = "vegetative"
            s.DOE = day
            # send signal to indicate crop emergence
            self._send_signal(signals.crop_emerged)
            
        elif s.STAGE == "vegetative":
            s.STAGE = "reproductive"
            s.DOA = day
                        
        elif s.STAGE == "reproductive":
            s.STAGE = "mature"
            s.DOM = day
            if p.CROP_END_TYPE in ["maturity","earliest"]:
                self._send_signal(signal=signals.crop_finish,
                                  day=day, finish="maturity",
                                  crop_delete=True)
        elif s.STAGE == "mature":
            msg = "Cannot move to next phenology stage: maturity already reached!"
            raise exc.PCSEError(msg)

        else: # Problem no stage defined
            msg = "No STAGE defined in phenology submodule."
            raise exc.PCSEError(msg)
        
        msg = "Changed phenological stage '%s' to '%s' on %s"
        self.logger.info(msg % (current_STAGE, s.STAGE, day))

    #---------------------------------------------------------------------------
    def _on_CROP_FINISH(self, day, finish_type=None):
        """Handler for setting day of harvest (DOH). Although DOH is not
        strictly related to phenology (but to management) this is the most
        logical place to put it.
        """
        if finish_type == 'harvest':
            self._for_finalize["DOH"] = day


class DVS_Phenology_Wrapper(SimulationObject):
    """This is wrapper class that wraps the DVS_phenology simulation object
    in such a way that it can run directly in the Engine. This means that
    it can be used direct in a configuration file::

        CROP = DVS_Phenology_Wrapper

    This is useful for running only the phenology instead of the entire
    crop simulation.
    """
    phenology = Instance(SimulationObject)

    def initialize(self, day, kiosk, cropdata, soildata, sitedata, start_type, stop_type):
        self.phenology = DVS_Phenology(day, kiosk, cropdata, start_type,  stop_type)

    def calc_rates(self, day, drv):
        self.phenology.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        self.phenology.integrate(day, delt)
