from pcse.base_classes import SimulationObject, ParamTemplate, StatesTemplate
from pcse.lintul.stateVariables import StateVariables
from pcse.lintul import lintul3lib
from pcse.lintul import lintul3
from pcse.lintul.lintul3lib import notNull, INSW, REAAND
from pcse.decorators import prepare_rates, prepare_states
from pcse.traitlets import Float
from numpy.ma.core import exp
from pcse.lintul.lintul3 import SubModel


class Lintul3Soil(SubModel):
    """
* ORIGINAL COPYRGIGHT NOTICE:    
*-------------------------------------------------------------------------*
* Copyright 2013. Wageningen University, Plant Production Systems group,  *
* P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
* You may not use this work except in compliance with the Licence.        *
* You may obtain a copy of the Licence at:                                *
*                                                                         *
* http://models.pps.wur.nl/content/licence-agreement                      *
*                                                                         *
* Unless required by applicable law or agreed to in writing, software     *
* distributed under the Licence is distributed on an "AS IS" basis,       *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
*-------------------------------------------------------------------------*

************************************************************************
*   LINTUL3 is an extended version of LINTUL1 (the version of LINTUL   *
*          for optimal growth conditions) and LINTUL2 includes a simple*
*          water balance for studying effects of drought. LINTUL2-N    *
*          includes N-limitation on crop growth. The latter program is *
*          called LINTUL3.                                             *
*          Test version for spring wheat, using parameters for spring  *
*          wheat for Flevoland                                         *
************************************************************************
    """
    
    class Parameters(ParamTemplate):
        DRATE   = Float(-99)
        IRRIGF  = Float(-99)
        ROOTDM  = Float(-99)
        RRDMAX  = Float(-99)
        WCFC    = Float(-99)
        WCI     = Float(-99)
        WCST    = Float(-99)
        WCSUBS  = Float(-99)
        WCWP    = Float(-99)
        WMFAC   = Float(-99)
    
    
    class Lintul3SoilStates(StateVariables):
        WA      = Float(-99.)
        TRUNOF  = Float(-99.)
        TTRAN   = Float(-99.)
        TEVAP   = Float(-99.)
        TDRAIN  = Float(-99.)
        TRAIN   = Float(-99.)
        TEXPLO  = Float(-99.)
        TIRRIG  = Float(-99.)
        
             
    
    class InitialValues(object):
        
        def __init__(self, parameters):
            # Read initial states

            # Initial amount of water present in the rooted depth at the start of
            # the calculations, based on the initial water content (in mm).
            ROOTDI= 0.1
            self.WAI  = 1000. * ROOTDI * parameters.WCI
            
        
       
        
        
    def __init__(self, day, kiosk, *args, **kwargs):
        self.find_subroutines()
        super(Lintul3Soil, self).__init__(day, kiosk, *args, **kwargs)
        
    
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)

        # Read initial states
        init                = self.InitialValues(self.params)
        self.initialValues  = init
        initialStates       = self.Lintul3SoilStates.initialValues() 

        # Initial amount of water present in the rooted depth at the start of
        # the calculations, based on the initial water content (in mm).
        initialStates["WA"] = init.WAI
        
        # Initialize state variables
        self.states = self.Lintul3SoilStates(kiosk, publish=["WA"], **initialStates)
        self.states.initialize()
        
        
        
    def find_subroutines(self):
        self.drunir = lintul3lib.DRUNIR
        self.penman = lintul3lib.PENMAN
        
        
    @prepare_rates
    def calc_rates(self, day, drv):

        # dynamic calculations
        p = self.params
        s = self.states
        i = self.initialValues
        
        DELT = 1 # ???
        
        
        ROOTD = self.kiosk["ROOTD"]
        EMERG = self.kiosk["EMERG"]
        EVAP = self.kiosk["EVAP"]
        TRAN = self.kiosk["TRAN"]
        
        # Variables supplied by the weather system
        RAIN             = drv.RAIN * 10 # cm  --> mm CORRECTION FOR NON-STANDARD cm in CABO-WEATHER
                
        
        #  Water content in the rootzone
        WC  = 0.001* s.WA /notNull(ROOTD)
                      
        RROOTD = min(p.RRDMAX * INSW(WC - p.WCWP, 0., 1.) * EMERG,  p.ROOTDM - ROOTD)
        
        # Calling the subroutine for rates of drainage, runoff and irrigation.
        DRAIN, RUNOFF, IRRIG = self.drunir(RAIN, EVAP, TRAN, p.IRRIGF, p.DRATE, 
                                           DELT, s.WA, ROOTD, p.WCFC, p.WCST, p.WMFAC)
        
        
        #  Exploration of water in soil when roots grow downward.
        EXPLOR = 1000. * RROOTD * p.WCSUBS
                
        RWA = (RAIN+EXPLOR+IRRIG)-(RUNOFF+TRAN+EVAP+DRAIN)

        s.rWA     = RWA   
        s.rTEXPLO = EXPLOR
        
        s.rTEVAP  = EVAP  
        s.rTTRAN  = TRAN  
        s.rTRUNOF = RUNOFF
        s.rTIRRIG = IRRIG 
        s.rTRAIN  = RAIN  
        s.rTDRAIN = DRAIN
        
        TIME =  day.timetuple().tm_yday
        self.doOutput(self, TIME, locals().copy())

        WATBAL = (s.WA + (s.TRUNOF + s.TTRAN + s.TEVAP + s.TDRAIN)    # @UnusedVariable
                         - (i.WAI + s.TRAIN + s.TEXPLO + s.TIRRIG))         


        
    @prepare_states
    def integrate(self, day):
        states = self.states
        delta = 1.
        
        for s in states.listIntegratedStates():
            rate = getattr(states, 'r' + s)
            state = getattr(states, s)
            newvalue = state + delta * rate
            setattr(states, s, newvalue)
        

        
        
        
        