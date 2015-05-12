from pcse.base_classes import StatesTemplate, VariableKiosk
from pcse.traitlets import Float




class StateVariables(StatesTemplate):
        
    rates = {}
    __initialized = False

            
        
    def __setattr__(self, name, value):
        if self.rates.has_key(name):
            self.rates[name] = value
        elif not self.__initialized:
            object.__setattr__(self, name, value)
        else:
            super(StateVariables, self).__setattr__(name, value)
            
            
    def __getattr__(self, name):
        if self.rates.has_key(name):
            return self.rates[name]
        else:
            object.__getattribute__(self, name)


    def initialize(self):
        self.rates = {}
        self.__initialized = True
        
        for s in self.__class__.listIntegratedStates():
            self.rates['r' + s] = 0.0


    @classmethod
    def listIntegratedStates(cls):
        return sorted([a for a in cls.__dict__ if isinstance(getattr(cls, a), Float) and not a.startswith('_')])



    @classmethod
    def initialValues(cls):
        return dict((a, 0.0) for a in cls.__dict__ if isinstance(getattr(cls, a), Float) and not a.startswith('_'))




class States(StateVariables):
    TSUM  = Float(-99.)
    LAI   = Float(-99.)
    ANLV  = Float(-99.)
    ANST  = Float(-99.)
    ANRT  = Float(-99.)
    ANSO  = Float(-99.)
    NUPTT = Float(-99.)
    TNSOIL= Float(-99.)
    NLOSSL= Float(-99.)
    NLOSSR= Float(-99.)
    WLVG  = Float(-99.)
    WLVD  = Float(-99.)
    WST   = Float(-99.)
    WSO   = Float(-99.)
    WRT   = Float(-99.)
    ROOTD = Float(-99.)
    GTSUM = Float(-99.)
    WDRT  = Float(-99.)
    CUMPAR= Float(-99.)
    WA    = Float(-99.)
    TEXPLO= Float(-99.)
    TEVAP = Float(-99.)
    TTRAN = Float(-99.)
    TRUNOF= Float(-99.)
    TIRRIG= Float(-99.)
    TRAIN = Float(-99.)
    TDRAIN= Float(-99.)
    
    

if (__name__ == "__main__"):
    
    kiosk  = VariableKiosk()
    init = States.initialValues()
    states = States(kiosk, publish=[], **init) 
    
    
    
    states.initialize()
    print states.listIntegratedStates()
    print states.rTSUM
          
    states.rTSUM = 123.
    
    print states.rTSUM
