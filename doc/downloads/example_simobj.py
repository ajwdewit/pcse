from .base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .decorators import prepare_rates, prepare_states
from .traitlets import Float

class CropProcess(SimulationObject):
    
    helper_variable1 = Float(-99.)
    
    class Parameters(ParamTemplate):
        PAR1  = Float(-99.)
        PAR2  = Float(-99.)
    
    class StateVariables(StatesTemplate):
        STATE1 = Float(-99.)
        STATE2 = Float(-99.)
        RATIO  = Float(-99.)

    class RateVariables(RatesTemplate):
        RATE1 = Float(-99.)
        RATE2 = Float(-99.)

    def initialize(self, day, kiosk, parametervalues):
        """Initializes the SimulationObject with given parametervalues."""

        self.kiosk = kiosk
        
        # Initialize parameters
        self.params = self.Parameters(parametervalues)

        # Initialize state variables
        self.states = self.StateVariables(kiosk, publish="STATE1", STATE1=0.,
                                          STATE2=10.)
        
        # Initialize rate variables
        self.rates = self.RateVariables(kiosk, publish=["RATE1","RATE2"])
        
    @prepare_rates
    def calc_rates(self, day, drv):
        """Calculate the rates of change."""
        
        self.rates.RATE1 = self.params.PAR1 * drv.IRRAD
        self.rates.RATE2 = self.params.PAR2 * (drv.TEMP)**2
        
        msg = "Rate calculation finished on CropProcess on %s." % day
        self.logger.info(msg)
    
    @prepare_states
    def integrate(self, day, delt):
        """Integrate the rates of change on the current state variables
        multiplied by the time-step
        """
        
        self.states.STATE1 += self.rates.RATE1 * delt
        self.states.STATE2 += self.rates.RATE2 * delts

        msg = "State update finished on CropProcess on %s." % day
        self.logger.info(msg)

    @prepare_states
    def finalize(self, day):
        try:
            self.states.RATIO = self.states.STATE1/self.states.STATE2
        except ZeroDivisionError:
            self.states.RATIO = -99.
            msg = "Failed to calculate STATE1/STATE2 ratio: division by zero"
            self.logger.warn(msg)

        # Ensure that embedded SimulationObjects get finalized as well
        SimulationObject.finalize(self, day)
