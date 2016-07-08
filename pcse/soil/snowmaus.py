# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from ..traitlets import Float
from ..decorators import prepare_rates, prepare_states
from ..util import limit, merge_dict
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals
from .. import exceptions as exc

class SnowMAUS(SimulationObject):
    """Simple snow accumulation model for agrometeorological applications.
    
    This is an implementation of the SnowMAUS model which describes the
    accumulation and melt of snow due to precipitation, snowmelt and
    sublimation. The SnowMAUS model is designed to keep track of the thickness
    of the layer of water that is present as snow on the surface, e.g. the
    Snow Water Equivalent Depth (state variable SWEDEPTH [cm]). Conversion of 
    the SWEDEPTH to the actual snow depth (state variable SNOWDEPTH [cm])
    is done by dividing the SWEDEPTH with the snow density in [cm_water/cm_snow].
    
    Snow density is taken as a fixed value despite the fact that the snow
    density is known to vary with the type of snowfall, the temperature and
    the age of the snow pack. However, more complicated algorithms for snow
    density would not be consistent with the simplicy of SnowMAUS.
    
    A drawback of the current implementation is that there is no link to the
    water balance yet.
    
    Reference:
    M. Trnka, E. Kocmánková, J. Balek, J. Eitzinger, F. Ruget, H. Formayer,
    P. Hlavinka, A. Schaumberger, V. Horáková, M. Možný, Z. Žalud,
    Simple snow cover model for agrometeorological applications,
    Agricultural and Forest Meteorology, Volume 150, Issues 7–8, 15 July 2010,
    Pages 1115-1127, ISSN 0168-1923

    http://dx.doi.org/10.1016/j.agrformet.2010.04.012

    **Simulation parameters:** (provide in crop, soil and sitedata dictionary)
    
    ============ =========================================== =======  ==========
     Name         Description                                 Type     Unit
    ============ =========================================== =======  ==========
    TMINACCU1    Upper critical minimum temperature for snow  SSi      |C|
                 accumulation. 
    TMINACCU2    Lower critical minimum temperature for snow  SSi      |C|
                 accumulation
    TMINCRIT     Critical minimum temperature for snow melt   SSi      |C|
    TMAXCRIT     Critical maximum temperature for snow melt   SSi      |C|
    RMELT        Melting rate per day per degree Celcius      SSi      |cmC-1day-1|
                 above the critical minimum temperature.
    SCTHRESHOLD  Snow water equivalent above which the        SSi      cm
                 sublimation is taken into account. 
    SNOWDENSITY  Density of snow                              SSi      cm/cm
    SWEDEPTHI    Initial depth of layer of water present as
                 snow on the soil surface                     SSi      cm
    ============ =========================================== =======  ==========
      
    **State variables:**

    =============== ========================================== ==== ============
     Name           Description                                Pbl      Unit
    =============== ========================================== ==== ============
    SWEDEPTH        Depth of layer of water present as snow     N    cm
                    on the surface
    SNOWDEPTH       Depth of snow present on the surface.       Y    cm
    =============== ========================================== ==== ============
    
    **Rate variables:**

    ============ ============================================= ==== ============
     Name        Description                                   Pbl     Unit
    ============ ============================================= ==== ============
    RSNOWACCUM   Rate of snow accumulation                      N    |cmday-1|
    RSNOWSUBLIM  Rate of snow sublimation                       N    |cmday-1|
    RSNOWMELT    Rate of snow melting                           N    |cmday-1|
    ============ ============================================= ==== ============
    """
    
    class Parameters(ParamTemplate):
        TMINACCU1 = Float(-99.)
        TMINACCU2 = Float(-99.)
        TMINCRIT  = Float(-99.)
        TMAXCRIT  = Float(-99.)
        RMELT     = Float(-99.)
        SCTHRESHOLD = Float(-99.)
        SNOWDENSITY = Float(-99)
        SWEDEPTHI   = Float(-99)


    class StateVariables(StatesTemplate):
        SWEDEPTH = Float(-99.)
        SNOWDEPTH = Float(-99.)

    class RateVariables(RatesTemplate):
        RSNOWACCUM  = Float(-99.)
        RSNOWSUBLIM = Float(-99.)
        RSNOWMELT   = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param sitedata: dictionary with WOFOST sitedata key/value pairs
        """
      
        if parvalues["SNOWDENSITY"] <= 0.:
            msg = ("SNOWDENSITY parameter of SnowMAUS module cannot <= zero " +
                   "to avoid division by zero")
            raise exc.ParameterError(msg)
        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk)

        SWEDEPTH = self.params.SWEDEPTHI
        SNOWDEPTH = SWEDEPTH/self.params.SNOWDENSITY
        self.states = self.StateVariables(kiosk, SWEDEPTH=SWEDEPTH, 
                                          SNOWDEPTH=SNOWDEPTH, publish="SNOWDEPTH")

    @prepare_rates
    def calc_rates(self, day, drv):
        p = self.params
        r = self.rates
        s = self.states
        
        # Snow accumulation rate (RSNOWACCUM)
        if drv.TMIN <= p.TMINACCU2:
            r.RSNOWACCUM = drv.RAIN
        elif drv.TMIN >= p.TMINACCU1:
            r.RSNOWACCUM = 0.
        else:
            rr = (drv.TMIN - p.TMINACCU2)/ abs(p.TMINACCU1 - p.TMINACCU2)
            r.RSNOWACCUM = (1-rr) * drv.RAIN

        # Snow sublimation rate (RSNOWSUBLIM)
        if s.SWEDEPTH > p.SCTHRESHOLD and r.RSNOWACCUM == 0.:
            RSNOWSUBLIM = drv.E0
        else:
            RSNOWSUBLIM = 0.
        # Avoid sublimating more snow than available
        r.RSNOWSUBLIM = limit(0, s.SWEDEPTH, RSNOWSUBLIM)        

        # Snow melting rate (RSNOWMELT)
        if drv.TMIN < p.TMINCRIT:
            RSNOWMELT = 0.
        else:
            if drv.TMIN <= 0. and drv.TMAX < p.TMAXCRIT:
                RSNOWMELT = 0.
            else:
                RSNOWMELT = (drv.TMIN - p.TMINCRIT) * p.RMELT

        # Avoid melting more snow than available
        r.RSNOWMELT = limit(0, (s.SWEDEPTH - r.RSNOWSUBLIM), RSNOWMELT)        

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        r = self.rates
        p = self.params
        
        s.SWEDEPTH += (r.RSNOWACCUM - r.RSNOWSUBLIM - r.RSNOWMELT)
        s.SNOWDEPTH = s.SWEDEPTH/p.SNOWDENSITY
