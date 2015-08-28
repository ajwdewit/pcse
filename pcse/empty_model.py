from pcse.base_classes import SimulationObject, ParamTemplate, AncillaryObject
from pcse.decorators import prepare_rates, prepare_states


class EmptyModel(SimulationObject):
    """
    dummy model: NoOp
    """


    def initialize(self, day, kiosk, *args, **kwargs):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        pass
        
        
    @prepare_rates
    def calc_rates(self, day, drv):
        pass
    
    
    
    @prepare_states
    def integrate(self, day):
        pass


        
class EmptyAncillaryObject(AncillaryObject):
    
    def initialize(self, day, kiosk, *args, **kwargs):
        pass
    
    def __call__ (self, *args, **kwargs):
        pass
