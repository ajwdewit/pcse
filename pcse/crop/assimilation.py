# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""SimulationObjects implementing |CO2| Assimilation for use with PCSE.
"""
from __future__ import print_function
from math import sqrt, exp, cos, pi

from ..traitlets import Instance, Float, AfgenTrait

from ..util import limit, astro, doy
from ..base_classes import ParamTemplate, SimulationObject

try:
    from ..futil import totass as ftotass
    from ..futil import astro as fastro
except ImportError as exc:
    ftotass = fastro = None


def totass(DAYL, AMAX, EFF, LAI, KDIF, AVRAD, DIFPP, DSINBE, SINLD, COSLD):
    """ This routine calculates the daily total gross CO2 assimilation by
    performing a Gaussian integration over time. At three different times of
    the day, irradiance is computed and used to calculate the instantaneous
    canopy assimilation, whereafter integration takes place. More information
    on this routine is given by Spitters et al. (1988).

    FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
    name   type meaning                                    units  class
    ----   ---- -------                                    -----  -----
    DAYL    R4  Astronomical daylength (base = 0 degrees)     h      I
    AMAX    R4  Assimilation rate at light saturation      kg CO2/   I
                                                          ha leaf/h
    EFF     R4  Initial light use efficiency              kg CO2/J/  I
                                                          ha/h m2 s
    LAI     R4  Leaf area index                             ha/ha    I
    KDIF    R4  Extinction coefficient for diffuse light             I
    AVRAD   R4  Daily shortwave radiation                  J m-2 d-1 I
    DIFPP   R4  Diffuse irradiation perpendicular to direction of
                light                                      J m-2 s-1 I
    DSINBE  R4  Daily total of effective solar height         s      I
    SINLD   R4  Seasonal offset of sine of solar height       -      I
    COSLD   R4  Amplitude of sine of solar height             -      I
    DTGA    R4  Daily total gross assimilation           kg CO2/ha/d O

    Authors: Daniel van Kraalingen
    Date   : April 1991

    Python version:
    Authors: Allard de Wit
    Date   : September 2011
    """

    # Gauss points and weights
    XGAUSS = [0.1127017, 0.5000000, 0.8872983]
    WGAUSS = [0.2777778, 0.4444444, 0.2777778]

    # calculation of assimilation is done only when it will not be zero
    # (AMAX >0, LAI >0, DAYL >0)
    DTGA = 0.
    if (AMAX > 0. and LAI > 0. and DAYL > 0.):
        for i in range(3):
            HOUR   = 12.0+0.5*DAYL*XGAUSS[i]
            SINB   = max(0.,SINLD+COSLD*cos(2.*pi*(HOUR+12.)/24.))
            PAR    = 0.5*AVRAD*SINB*(1.+0.4*SINB)/DSINBE
            PARDIF = min(PAR,SINB*DIFPP)
            PARDIR = PAR-PARDIF
            FGROS = assim(AMAX,EFF,LAI,KDIF,SINB,PARDIR,PARDIF)
            DTGA += FGROS*WGAUSS[i]
    DTGA *= DAYL

    return DTGA


def assim(AMAX, EFF, LAI, KDIF, SINB, PARDIR, PARDIF):
    """This routine calculates the gross CO2 assimilation rate of
    the whole crop, FGROS, by performing a Gaussian integration
    over depth in the crop canopy. At three different depths in
    the canopy, i.e. for different values of LAI, the
    assimilation rate is computed for given fluxes of photosynthe-
    tically active radiation, whereafter integration over depth
    takes place. More information on this routine is given by
    Spitters et al. (1988). The input variables SINB, PARDIR
    and PARDIF are calculated in routine TOTASS.

    Subroutines and functions called: none.
    Called by routine TOTASS.

    Author: D.W.G. van Kraalingen, 1986

    Python version:
    Allard de Wit, 2011
    """
    # Gauss points and weights
    XGAUSS = [0.1127017, 0.5000000, 0.8872983]
    WGAUSS = [0.2777778, 0.4444444, 0.2777778]

    SCV = 0.2

    # 13.2 extinction coefficients KDIF, KDIRBL, KDIRT
    REFH = (1.-sqrt(1.-SCV))/(1.+sqrt(1.-SCV))
    REFS = REFH*2./(1.+1.6*SINB)
    KDIRBL = (0.5/SINB)*KDIF/(0.8*sqrt(1.-SCV))
    KDIRT = KDIRBL*sqrt(1.-SCV)

    #13.3 three-point Gaussian integration over LAI
    FGROS = 0.
    for i in range(3):
        LAIC = LAI*XGAUSS[i]
        # absorbed diffuse radiation (VISDF),light from direct
        # origine (VIST) and direct light (VISD)
        VISDF  = (1.-REFS)*PARDIF*KDIF  *exp(-KDIF  *LAIC)
        VIST   = (1.-REFS)*PARDIR*KDIRT *exp(-KDIRT *LAIC)
        VISD   = (1.-SCV) *PARDIR*KDIRBL*exp(-KDIRBL*LAIC)

        # absorbed flux in W/m2 for shaded leaves and assimilation
        VISSHD = VISDF+VIST-VISD
        FGRSH  = AMAX*(1.-exp(-VISSHD*EFF/max(2.0,AMAX)))

        # direct light absorbed by leaves perpendicular on direct
        # beam and assimilation of sunlit leaf area
        VISPP  = (1.-SCV)*PARDIR/SINB
        if (VISPP <= 0.):
            FGRSUN = FGRSH
        else:
            FGRSUN = AMAX*(1.-(AMAX-FGRSH) \
                     *(1.-exp(-VISPP*EFF/max(2.0,AMAX)))/ (EFF*VISPP))

        # fraction of sunlit leaf area (FSLLA) and local
        # assimilation rate (FGL)
        FSLLA  = exp(-KDIRBL*LAIC)
        FGL    = FSLLA*FGRSUN+(1.-FSLLA)*FGRSH

        # integration
        FGROS += FGL*WGAUSS[i]

    FGROS  = FGROS*LAI
    return FGROS


class WOFOST_Assimilation(SimulationObject):
    """Class implementing a WOFOST/SUCROS style assimilation routine.
    
    WOFOST calculates the daily gross |CO2| assimilation rate of a crop
    from the absorbed radiation and the photosynthesis-light response curve
    of individual leaves. This response is dependent on temperature and
    leaf age. The absorbed radiation is calculated from the total incoming
    radiation and the leaf area. Daily gross |CO2| assimilation is obtained
    by integrating the assimilation rates over the leaf layers and over the
    day.
      
    *Simulation parameters*
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    AMAXTB   Max. leaf |CO2| assim. rate as a function of   TCr     |kg ha-1 hr-1|
             of DVS
    EFFTB    Light use effic. single leaf as a function     TCr     |kg ha-1 hr-1 /(J m-2 s-1)|
             of daily mean temperature                                    
    KDIFTB   Extinction coefficient for diffuse visible     TCr      -
             as function of DVS
    TMPFTB   Reduction factor of AMAX as function of        TCr      -
             daily mean temperature.
    TMNFTB   Reduction factor of AMAX as function of        TCr      -
             daily minimum temperature.
    =======  ============================================= =======  ============
    
    *State and rate variables*
    
    `WOFOST_Assimilation` has no state/rate variables, but calculates the
    rate of assimilation which is returned directly from the `__call__()`
    method.
    
    *Signals sent or handled*
    
    None
        
    
    *External dependencies:*
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    LAI      Leaf area index                     Leaf_dynamics       -
    =======  =================================== =================  ============
    """
    
    class Parameters(ParamTemplate):
        AMAXTB = AfgenTrait()
        EFFTB  = AfgenTrait()
        KDIFTB = AfgenTrait()
        TMPFTB = AfgenTrait()
        TMNFTB = AfgenTrait()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        :returns: the assimilation rate in |kg ha-1 d-1| using __call__()
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk
    
    def __call__(self, day, drv):
        # Check if fortran versions can be used otherwise use native python
        if ftotass is not None:
            PGASS = self.___call__fortran(day, drv)
        else:
            PGASS = self.__call__python(day, drv)

        return PGASS

    def ___call__fortran(self, day, drv):
        """Calls fortran versions of ASTRO and TOTASS
        """
        params = self.params

        # published states from the kiosk
        DVS = self.kiosk["DVS"]
        LAI = self.kiosk["LAI"]

        # 2.19  photoperiodic daylength
        IDAY = doy(day)
        DAYL, DAYLP, SINLD, COSLD, DIFPP, ATMTR, DSINBE = \
             fastro(IDAY, drv.LAT, drv.IRRAD)

        # 2.20  daily dry matter production

        # gross assimilation and correction for sub-optimum
        # average day temperature
        AMAX = params.AMAXTB(DVS)
        AMAX *= params.TMPFTB(drv.DTEMP)
        KDIF = params.KDIFTB(DVS)
        EFF  = params.EFFTB(drv.DTEMP)
        DTGA = ftotass(DAYL, AMAX, EFF, LAI, KDIF, drv.IRRAD, DIFPP,
                       DSINBE, SINLD, COSLD)

        # correction for low minimum temperature potential
        # assimilation in kg CH2O per ha
        DTGA *= params.TMNFTB(drv.TMINRA)
        PGASS = DTGA * 30./44.

        return PGASS

    def __call__python(self, day, drv):
        params = self.params

        # published states from the kiosk
        DVS = self.kiosk["DVS"]
        LAI = self.kiosk["LAI"]
        
        # 2.19  photoperiodic daylength
        DAYL,DAYLP,SINLD,COSLD,DIFPP,ATMTR,DSINBE,ANGOT = astro(day, drv.LAT, drv.IRRAD)

        # 2.20  daily dry matter production

        # gross assimilation and correction for sub-optimum
        # average day temperature
        AMAX = params.AMAXTB(DVS)
        AMAX *= params.TMPFTB(drv.DTEMP)
        KDIF = params.KDIFTB(DVS)
        EFF  = params.EFFTB(drv.DTEMP)        
        DTGA = totass(DAYL, AMAX, EFF, LAI, KDIF, drv.IRRAD, DIFPP,
                      DSINBE, SINLD, COSLD)

        # correction for low minimum temperature potential
        # assimilation in kg CH2O per ha
        DTGA *= params.TMNFTB(drv.TMINRA)
        PGASS = DTGA * 30./44.
        
        return PGASS


class WOFOST_Assimilation2(SimulationObject):
    """Class implementing a WOFOST/SUCROS style assimilation routine including
    effect of changes in atmospheric CO2 concentration.

    WOFOST calculates the daily gross |CO2| assimilation rate of a crop
    from the absorbed radiation and the photosynthesis-light response curve
    of individual leaves. This response is dependent on temperature and
    leaf age. The absorbed radiation is calculated from the total incoming
    radiation and the leaf area. Daily gross |CO2| assimilation is obtained
    by integrating the assimilation rates over the leaf layers and over the
    day.



    *Simulation parameters* (To be provided in cropdata dictionary):

    =========  ============================================= =======  ============
     Name       Description                                   Type     Unit
    =========  ============================================= =======  ============
    AMAXTB     Max. leaf |CO2| assim. rate as a function of   TCr     |kg ha-1 hr-1|
               of DVS
    EFFTB      Light use effic. single leaf as a function     TCr     |kg ha-1 hr-1 /(J m-2 s-1)|
               of daily mean temperature
    KDIFTB     Extinction coefficient for diffuse visible     TCr      -
               as function of DVS
    TMPFTB     Reduction factor of AMAX as function of        TCr      -
               daily mean temperature.
    TMPFTB     Reduction factor of AMAX as function of        TCr      -
               daily minimum temperature.
    CO2AMAXTB  Correction factor for AMAX given atmos-        TCr      -
               pheric CO2 concentration.
    CO2EFFTB   Correction factor for EFF given atmos-         TCr      -
               pheric CO2 concentration.
    CO2        Atmopheric CO2 concentration                   SCr      ppm
    =========  ============================================= =======  ============

    *State and rate variables*

    `WOFOST_Assimilation2` has no state/rate variables, but calculates the
    rate of assimilation which is returned directly from the `__call__()`
    method.

    *Signals sent or handled*

    None


    *External dependencies:*

    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    LAI      Leaf area index                     Leaf_dynamics       -
    =======  =================================== =================  ============
    """

    class Parameters(ParamTemplate):
        AMAXTB = AfgenTrait()
        EFFTB = AfgenTrait()
        KDIFTB = AfgenTrait()
        TMPFTB = AfgenTrait()
        TMNFTB = AfgenTrait()
        CO2AMAXTB = AfgenTrait()
        CO2EFFTB = AfgenTrait()
        CO2 = Float(-99.)

    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this Engine instance
        :param cropdata: dictionary with cropdata key/value pairs
        :returns: the assimilation rate using __call__()
        """

        self.params = self.Parameters(cropdata)
        self.kiosk  = kiosk

    def __call__(self, day, drv):
        params = self.params

        # published states from the kiosk
        DVS = self.kiosk["DVS"]
        LAI = self.kiosk["LAI"]

        # 2.19  photoperiodic daylength
        DAYL,DAYLP,SINLD,COSLD,DIFPP,ATMTR,DSINBE,ANGOT = astro(day, drv.LAT, drv.IRRAD)

        # 2.20  daily dry matter production

        # gross assimilation and correction for sub-optimum
        # average day temperature
        AMAX = params.AMAXTB(DVS)

        # added IS june 2013
        # correction for changing CO2 concentration
        AMAX *= params.CO2AMAXTB(params.CO2)
        AMAX *= params.TMPFTB(drv.DTEMP)
        KDIF = params.KDIFTB(DVS)
        EFF  = params.EFFTB(drv.DTEMP)*params.CO2EFFTB(params.CO2)
        DTGA = totass(DAYL, AMAX, EFF, LAI, KDIF, drv.IRRAD, DIFPP,
                      DSINBE, SINLD, COSLD)

        # correction for low minimum temperature potential
        # assimilation in kg CH2O per ha
        DTGA *= params.TMNFTB(drv.TMINRA)
        PGASS = DTGA * 30./44.

        return PGASS

