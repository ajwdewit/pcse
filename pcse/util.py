# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Miscellaneous utilities for PCSE
"""
import os, sys
import datetime
import copy
from math import log10, cos, sin, asin, sqrt, exp
from collections import namedtuple
from bisect import bisect_left
try:
    # support both python2 and python3
    from collections import MutableMapping
except ImportError:
    from UserDict import DictMixin as MutableMapping
import textwrap
import sqlite3
import pdb


import numpy as np

from . import exceptions as exc

Celsius2Kelvin = lambda x: x + 273.16
hPa2kPa = lambda x: x/10.
# Saturated Vapour pressure [kPa] at temperature temp [C]
SatVapourPressure = lambda temp: 0.6108 * exp((17.27 * temp) / (237.3 + temp))

def reference_ET(DAY, LAT, ELEV, TMIN, TMAX, IRRAD, VAP, WIND,
                 ANGSTA, ANGSTB, ETMODEL="PM", **kwargs):
    """Calculates reference evapotranspiration values E0, ES0 and ET0.

    The open water (E0) and bare soil evapotranspiration (ES0) are calculated with
    the modified Penman approach, while the references canopy evapotranspiration is
    calculated with the modified Penman or the Penman-Monteith approach, the latter
    is the default.

    Input variables::

        DAY     -  Python datetime.date object                      -
        LAT     -  Latitude of the site                          degrees
        ELEV    -  Elevation above sea level                        m
        TMIN    -  Minimum temperature                              C
        TMAX    -  Maximum temperature                              C
        IRRAD   -  Daily shortwave radiation                     J m-2 d-1
        VAP     -  24 hour average vapour pressure                 hPa
        WIND    -  24 hour average windspeed at 2 meter            m/s
        ANGSTA  -  Empirical constant in Angstrom formula           -
        ANGSTB  -  Empirical constant in Angstrom formula           -
        ETMODEL -  Indicates if the canopy reference ET should     PM|P
                   be calculated with the Penman-Monteith method
                   (PM) or the modified Penman method (P)

    Output is a tuple (E0, ES0, ET0)::

        E0      -  Penman potential evaporation from a free
                   water surface [mm/d]
        ES0     -  Penman potential evaporation from a moist
                   bare soil surface [mm/d]
        ET0     -  Penman or Penman-Monteith potential evapotranspiration from a
                   crop canopy [mm/d]

.. note:: The Penman-Monteith algorithm is valid only for a reference canopy and
    therefore it is not used to calculate the reference values for bare soil and
    open water (ES0, E0).

    The background is that the Penman-Monteith model is basically a surface
    energy balance where the net solar radiation is partitioned over latent and
    sensible heat fluxes (ignoring the soil heat flux). To estimate this
    partitioning, the PM method makes a connection between the surface
    temperature and the air temperature. However, the assumptions
    underlying the PM model are valid only when the surface where this
    partitioning takes place is the same for the latent and sensible heat
    fluxes.

    For a crop canopy this assumption is valid because the leaves of the
    canopy form the surface where both latent heat flux (through stomata)
    and sensible heat flux (through leaf temperature) are partitioned.
    For a soil, this principle does not work because when the soil is
    drying the evaporation front will quickly disappear below the surface
    and therefore the assumption that the partitioning surface is the
    same does not hold anymore.

    For water surfaces, the assumptions underlying PM do not hold
    because there is no direct relationship between the temperature
    of the water surface and the net incoming radiation as radiation is
    absorbed by the water column and the temperature of the water surface
    is co-determined by other factors (mixing, etc.). Only for a very
    shallow layer of water (1 cm) the PM methodology could be applied.

    For bare soil and open water the Penman model is preferred. Although it
    partially suffers from the same problems, it is calibrated somewhat
    better for open water and bare soil based on its empirical wind
    function.

    Finally, in crop simulation models the open water evaporation and
    bare soil evaporation only play a minor role (pre-sowing conditions
    and flooded rice at early stages), it is not worth investing much
    effort in improved estimates of reference value for E0 and ES0.
    """
    if ETMODEL not in ["PM", "P"]:
        msg = "Variable ETMODEL can have values 'PM'|'P' only."
        raise RuntimeError(msg)

    E0, ES0, ET0 = penman(DAY, LAT, ELEV, TMIN, TMAX, IRRAD, VAP, WIND,
                          ANGSTA, ANGSTB)
    if ETMODEL == "PM":
        ET0 = penman_monteith(DAY, LAT, ELEV, TMIN, TMAX, IRRAD, VAP, WIND)

    return E0, ES0, ET0


def penman(DAY, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND2, ANGSTA, ANGSTB):
    """Calculates E0, ES0, ET0 based on the Penman model.
    
     This routine calculates the potential evapo(transpi)ration rates from
     a free water surface (E0), a bare soil surface (ES0), and a crop canopy
     (ET0) in mm/d. For these calculations the analysis by Penman is followed
     (Frere and Popov, 1979;Penman, 1948, 1956, and 1963).
     Subroutines and functions called: ASTRO, LIMIT.

    Input variables::
    
        DAY     -  Python datetime.date object                                    -
        LAT     -  Latitude of the site                        degrees   
        ELEV    -  Elevation above sea level                      m      
        TMIN    -  Minimum temperature                            C
        TMAX    -  Maximum temperature                            C      
        AVRAD   -  Daily shortwave radiation                   J m-2 d-1 
        VAP     -  24 hour average vapour pressure               hPa     
        WIND2   -  24 hour average windspeed at 2 meter          m/s     
        ANGSTA  -  Empirical constant in Angstrom formula         -
        ANGSTB  -  Empirical constant in Angstrom formula         -

    Output is a tuple (E0,ES0,ET0)::
    
        E0      -  Penman potential evaporation from a free water surface [mm/d]
        ES0     -  Penman potential evaporation from a moist bare soil surface [mm/d]
        ET0     -  Penman potential transpiration from a crop canopy [mm/d]
    """
    # psychrometric instrument constant (mbar/Celsius-1)
    # albedo for water surface, soil surface and canopy
    # latent heat of evaporation of water (J/kg=J/mm)
    # Stefan Boltzmann constant (in J/m2/d/K4, e.g multiplied by 24*60*60)
    PSYCON = 0.67; REFCFW = 0.05; REFCFS = 0.15; REFCFC = 0.25
    LHVAP = 2.45E6; STBC =  5.670373E-8 * 24*60*60 # (=4.9E-3)

    # preparatory calculations
    # mean daily temperature and temperature difference (Celsius)
    # coefficient Bu in wind function, dependent on temperature
    # difference
    TMPA = (TMIN+TMAX)/2.
    TDIF = TMAX - TMIN
    BU = 0.54 + 0.35 * limit(0.,1.,(TDIF-12.)/4.)

    # barometric pressure (mbar)
    # psychrometric constant (mbar/Celsius)
    PBAR = 1013.*exp(-0.034*ELEV/(TMPA+273.))
    GAMMA = PSYCON*PBAR/1013.

    # saturated vapour pressure according to equation of Goudriaan
    # (1977) derivative of SVAP with respect to temperature, i.e. 
    # slope of the SVAP-temperature curve (mbar/Celsius);
    # measured vapour pressure not to exceed saturated vapour pressure

    SVAP = 6.10588 * exp(17.32491*TMPA/(TMPA+238.102))
    DELTA = 238.102*17.32491*SVAP/(TMPA+238.102)**2
    VAP = min(VAP, SVAP)

    # the expression n/N (RELSSD) from the Penman formula is estimated
    # from the Angstrom formula: RI=RA(A+B.n/N) -> n/N=(RI/RA-A)/B,
    # where RI/RA is the atmospheric transmission obtained by a CALL
    # to ASTRO:

    r = astro(DAY, LAT, AVRAD)
    RELSSD = limit(0., 1., (r.ATMTR-abs(ANGSTA))/abs(ANGSTB))

    # Terms in Penman formula, for water, soil and canopy

    # net outgoing long-wave radiation (J/m2/d) acc. to Brunt (1932)
    RB = STBC*(TMPA+273.)**4*(0.56-0.079*sqrt(VAP))*(0.1+0.9*RELSSD)

    # net absorbed radiation, expressed in mm/d
    RNW = (AVRAD*(1.-REFCFW)-RB)/LHVAP
    RNS = (AVRAD*(1.-REFCFS)-RB)/LHVAP
    RNC = (AVRAD*(1.-REFCFC)-RB)/LHVAP

    # evaporative demand of the atmosphere (mm/d)
    EA  = 0.26 * max(0.,(SVAP-VAP)) * (0.5+BU*WIND2)
    EAC = 0.26 * max(0.,(SVAP-VAP)) * (1.0+BU*WIND2)

    # Penman formula (1948)
    E0  = (DELTA*RNW+GAMMA*EA)/(DELTA+GAMMA)
    ES0 = (DELTA*RNS+GAMMA*EA)/(DELTA+GAMMA) 
    ET0 = (DELTA*RNC+GAMMA*EAC)/(DELTA+GAMMA)

    # Ensure reference evaporation >= 0.
    E0  = max(0., E0)
    ES0 = max(0., ES0)
    ET0 = max(0., ET0)
    
    return E0, ES0, ET0


def penman_monteith(DAY, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND2):
    """Calculates reference ET0 based on the Penman-Monteith model.

     This routine calculates the potential evapotranspiration rate from
     a reference crop canopy (ET0) in mm/d. For these calculations the
     analysis by FAO is followed as laid down in the FAO publication
     `Guidelines for computing crop water requirements - FAO Irrigation
     and drainage paper 56 <http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents>`_

    Input variables::

        DAY   -  Python datetime.date object                   -
        LAT   -  Latitude of the site                        degrees
        ELEV  - Elevation above sea level                      m
        TMIN  - Minimum temperature                            C
        TMAX  - Maximum temperature                            C
        AVRAD - Daily shortwave radiation                   J m-2 d-1
        VAP   - 24 hour average vapour pressure               hPa
        WIND2 - 24 hour average windspeed at 2 meter          m/s

    Output is:

        ET0   - Penman-Monteith potential transpiration
                rate from a crop canopy                     [mm/d]
    """

    # psychrometric instrument constant (kPa/Celsius)
    PSYCON = 0.665
    # albedo and surface resistance [sec/m] for the reference crop canopy
    REFCFC = 0.23; CRES = 70.
    # latent heat of evaporation of water [J/kg == J/mm] and
    LHVAP = 2.45E6
    # Stefan Boltzmann constant (J/m2/d/K4, e.g multiplied by 24*60*60)
    STBC = 4.903E-3
    # Soil heat flux [J/m2/day] explicitly set to zero
    G = 0.

    # mean daily temperature (Celsius)
    TMPA = (TMIN+TMAX)/2.

    # Vapour pressure to kPa
    VAP = hPa2kPa(VAP)

    # atmospheric pressure (kPa)
    T = Celsius2Kelvin(TMPA)
    PATM = 101.3 * pow((T - (0.0065*ELEV))/T, 5.26)

    # psychrometric constant (kPa/Celsius)
    GAMMA = PSYCON * PATM * 1.0E-3

    # Derivative of SVAP with respect to mean temperature, i.e.
    # slope of the SVAP-temperature curve (kPa/Celsius);
    SVAP_TMPA = SatVapourPressure(TMPA)
    DELTA = (4098. * SVAP_TMPA)/pow((TMPA + 237.3), 2)

    # Daily average saturated vapour pressure [kPa] from min/max temperature
    SVAP_TMAX = SatVapourPressure(TMAX)
    SVAP_TMIN = SatVapourPressure(TMIN)
    SVAP = (SVAP_TMAX + SVAP_TMIN) / 2.

    # measured vapour pressure not to exceed saturated vapour pressure
    VAP = min(VAP, SVAP)

    # Longwave radiation according at Tmax, Tmin (J/m2/d)
    # and preliminary net outgoing long-wave radiation (J/m2/d)
    STB_TMAX = STBC * pow(Celsius2Kelvin(TMAX), 4)
    STB_TMIN = STBC * pow(Celsius2Kelvin(TMIN), 4)
    RNL_TMP = ((STB_TMAX + STB_TMIN) / 2.) * (0.34 - 0.14 * sqrt(VAP))

    # Clear Sky radiation [J/m2/DAY] from Angot TOA radiation
    # the latter is found through a call to astro()
    r = astro(DAY, LAT, AVRAD)
    CSKYRAD = (0.75 + (2e-05 * ELEV)) * r.ANGOT

    if CSKYRAD > 0:
        # Final net outgoing longwave radiation [J/m2/day]
        RNL = RNL_TMP * (1.35 * (AVRAD/CSKYRAD) - 0.35)

        # radiative evaporation equivalent for the reference surface
        # [mm/DAY]
        RN = ((1-REFCFC) * AVRAD - RNL)/LHVAP

        # aerodynamic evaporation equivalent [mm/day]
        EA = ((900./(TMPA + 273)) * WIND2 * (SVAP - VAP))

        # Modified psychometric constant (gamma*)[kPa/C]
        MGAMMA = GAMMA * (1. + (CRES/208.*WIND2))

        # Reference ET in mm/day
        ET0 = (DELTA * (RN-G))/(DELTA + MGAMMA) + (GAMMA * EA)/(DELTA + MGAMMA)
        ET0 = max(0., ET0)
    else:
        ET0 = 0.

    return ET0

#-------------------------------------------------------------------------------
def check_angstromAB(xA, xB):
    """Routine checks validity of Angstrom coefficients.
    
    This is the  python version of the FORTRAN routine 'WSCAB' in 'weather.for'.
    """
    MIN_A = 0.1
    MAX_A = 0.4
    MIN_B = 0.3
    MAX_B = 0.7
    MIN_SUM_AB = 0.6
    MAX_SUM_AB = 0.9
    A = abs(xA)
    B = abs(xB)
    SUM_AB = A + B
    if A < MIN_A or A > MAX_A:
        msg = "invalid Angstrom A value!"
        raise exc.PCSEError(msg)
    if B < MIN_B or B > MAX_B:
        msg = "invalid Angstrom B value!"
        raise exc.PCSEError(msg)
    if SUM_AB < MIN_SUM_AB or SUM_AB > MAX_SUM_AB:
        msg = "invalid sum of Angstrom A & B values!"
        raise exc.PCSEError(msg)

    return [A,B]


def wind10to2(wind10):
    """Converts windspeed at 10m to windspeed at 2m using log. wind profile
    """
    wind2 = wind10 * (log10(2./0.033) / log10(10/0.033))
    return wind2


def ea_from_tdew(tdew):
    """
    Calculates actual vapour pressure, ea [kPa] from the dewpoint temperature
    using equation (14) in the FAO paper. As the dewpoint temperature is the
    temperature to which air needs to be cooled to make it saturated, the
    actual vapour pressure is the saturation vapour pressure at the dewpoint
    temperature. This method is preferable to calculating vapour pressure from
    minimum temperature.

    Taken from fao_et0.py written by Mark Richards

    Reference:
    Allen, R.G., Pereira, L.S., Raes, D. and Smith, M. (1998) Crop
        evapotranspiration. Guidelines for computing crop water requirements,
        FAO irrigation and drainage paper 56)

    Arguments:
    tdew - dewpoint temperature [deg C]
    """
    # Raise exception:
    if (tdew < -95.0 or tdew > 65.0):
        # Are these reasonable bounds?
        msg = 'tdew=%g is not in range -95 to +60 deg C' % tdew
        raise ValueError(msg)

    tmp = (17.27 * tdew) / (tdew + 237.3)
    ea = 0.6108 * exp(tmp)
    return ea


def angstrom(day, latitude, ssd, cA, cB):
    """Compute global radiation using the Angstrom equation.
    
    Global radiation is derived from sunshine duration using the Angstrom
    equation:
    globrad = Angot * (cA + cB * (sunshine / daylength)
    
    :param day: day of observation (date object)
    :param latitude: Latitude of the observation
    :param ssd: Observed sunshine duration
    :param cA: Angstrom A parameter
    :param cB: Angstrom B parameter    
    :returns: the global radiation in J/m2/day
    """
    r = astro(day, latitude, 0)
    globrad = r.ANGOT * (cA + cB * (ssd / r.DAYL))
    return globrad
                       

def doy(day):
    """Converts a date or datetime object to day-of-year (Jan 1st = doy 1)
    """
    # Check if day is a date or datetime object
    if isinstance(day, datetime.date):
        pass
    elif isinstance(day, datetime.datetime):
        day = day.date()
    else:
        msg = "Parameter day is not a date or datetime object."
        raise RuntimeError(msg)

    return (day-datetime.date(day.year,1,1)).days + 1


def limit(min, max, v):
    """limits the range of v between min and max
    """

    if min > max:
        raise RuntimeError("Min value (%f) larger than max (%f)" % (min, max))
    
    if v < min:       # V below range: return min
        return min
    elif v < max:     # v within range: return v
        return v
    else:             # v above range: return max
        return max


def daylength(day, latitude, angle=-4, _cache={}):
    """Calculates the daylength for a given day, altitude and base.

    :param day:         date/datetime object
    :param latitude:    latitude of location
    :param angle:       The photoperiodic daylength starts/ends when the sun
        is `angle` degrees under the horizon. Default is -4 degrees.
    
    Derived from the WOFOST routine ASTRO.FOR and simplified to include only
    daylength calculation. Results are being cached for performance
    """
    #from unum.units import h

    # Check for range of latitude
    if abs(latitude) > 90.:
        msg = "Latitude not between -90 and 90"
        raise RuntimeError(msg)
    
    # Calculate day-of-year from date object day
    IDAY = doy(day)
    
    # Test if daylength for given (day, latitude, angle) was already calculated
    # in a previous run. If not (e.g. KeyError) calculate the daylength, store
    # in cache and return the value.
    try:
        return _cache[(IDAY, latitude, angle)]
    except KeyError:
        pass
    
    # constants
    RAD = 0.0174533
    PI = 3.1415926

    # map python functions to capitals
    SIN = sin
    COS = cos
    ASIN = asin
    REAL = float

    # calculate daylength
    ANGLE = angle
    LAT = latitude
    DEC = -ASIN(SIN(23.45*RAD)*COS(2.*PI*(REAL(IDAY)+10.)/365.))
    SINLD = SIN(RAD*LAT)*SIN(DEC)
    COSLD = COS(RAD*LAT)*COS(DEC)
    AOB   = (-SIN(ANGLE*RAD)+SINLD)/COSLD

    # daylength
    if abs(AOB) <= 1.0:
        DAYLP = 12.0*(1.+2.*ASIN((-SIN(ANGLE*RAD)+SINLD)/COSLD)/PI)
    elif AOB > 1.0:
        DAYLP = 24.0
    else:
        DAYLP =  0.0

    # store results in cache
    _cache[(IDAY, latitude, angle)] = DAYLP
    
    return DAYLP


def astro(day, latitude, radiation, _cache={}):
    """python version of ASTRO routine by Daniel van Kraalingen.
    
    This subroutine calculates astronomic daylength, diurnal radiation
    characteristics such as the atmospheric transmission, diffuse radiation etc.

    :param day:         date/datetime object
    :param latitude:    latitude of location
    :param radiation:   daily global incoming radiation (J/m2/day)

    output is a `namedtuple` in the following order and tags::

        DAYL      Astronomical daylength (base = 0 degrees)     h      
        DAYLP     Astronomical daylength (base =-4 degrees)     h      
        SINLD     Seasonal offset of sine of solar height       -      
        COSLD     Amplitude of sine of solar height             -      
        DIFPP     Diffuse irradiation perpendicular to
                  direction of light                         J m-2 s-1 
        ATMTR     Daily atmospheric transmission                -      
        DSINBE    Daily total of effective solar height         s
        ANGOT     Angot radiation at top of atmosphere       J m-2 d-1
 
    Authors: Daniel van Kraalingen
    Date   : April 1991
 
    Python version
    Author      : Allard de Wit
    Date        : January 2011
    """

    # Check for range of latitude
    if abs(latitude) > 90.:
        msg = "Latitude not between -90 and 90"
        raise RuntimeError(msg)
    LAT = latitude
        
    # Determine day-of-year (IDAY) from day
    IDAY = doy(day)
    
    # reassign radiation
    AVRAD = radiation

    # Test if variables for given (day, latitude, radiation) were already calculated
    # in a previous run. If not (e.g. KeyError) calculate the variables, store
    # in cache and return the value.
    try:
        return _cache[(IDAY, LAT, AVRAD)]
    except KeyError:
        pass

    # constants
    RAD = 0.0174533
    PI = 3.1415926
    ANGLE = -4.

    # map python functions to capitals
    SIN = sin
    COS = cos
    ASIN = asin
    REAL = float
    SQRT = sqrt
    ABS = abs

    # Declination and solar constant for this day
    DEC = -ASIN(SIN(23.45*RAD)*COS(2.*PI*(REAL(IDAY)+10.)/365.))
    SC  = 1370.*(1.+0.033*COS(2.*PI*REAL(IDAY)/365.))

    # calculation of daylength from intermediate variables
    # SINLD, COSLD and AOB
    SINLD = SIN(RAD*LAT)*SIN(DEC)
    COSLD = COS(RAD*LAT)*COS(DEC)
    AOB = SINLD/COSLD

    # For very high latitudes and days in summer and winter a limit is
    # inserted to avoid math errors when daylength reaches 24 hours in 
    # summer or 0 hours in winter.

    # Calculate solution for base=0 degrees
    if abs(AOB) <= 1.0:
        DAYL  = 12.0*(1.+2.*ASIN(AOB)/PI)
        # integrals of sine of solar height
        DSINB  = 3600.*(DAYL*SINLD+24.*COSLD*SQRT(1.-AOB**2)/PI)
        DSINBE = 3600.*(DAYL*(SINLD+0.4*(SINLD**2+COSLD**2*0.5))+
                 12.*COSLD*(2.+3.*0.4*SINLD)*SQRT(1.-AOB**2)/PI)
    else:
        if AOB >  1.0: DAYL = 24.0
        if AOB < -1.0: DAYL = 0.0
        # integrals of sine of solar height	
        DSINB = 3600.*(DAYL*SINLD)
        DSINBE = 3600.*(DAYL*(SINLD+0.4*(SINLD**2+COSLD**2*0.5)))

    # Calculate solution for base=-4 (ANGLE) degrees
    AOB_CORR = (-SIN(ANGLE*RAD)+SINLD)/COSLD
    if abs(AOB_CORR) <= 1.0:
        DAYLP = 12.0*(1.+2.*ASIN(AOB_CORR)/PI)
    elif AOB_CORR > 1.0:
        DAYLP = 24.0
    elif AOB_CORR < -1.0:
        DAYLP = 0.0

    # extraterrestrial radiation and atmospheric transmission
    ANGOT = SC*DSINB
    # Check for DAYL=0 as in that case the angot radiation is 0 as well
    if DAYL > 0.0:
        ATMTR = AVRAD/ANGOT
    else:
        ATMTR = 0.

    # estimate fraction diffuse irradiation
    if (ATMTR > 0.75):
        FRDIF = 0.23
    elif (ATMTR <= 0.75) and (ATMTR > 0.35):
        FRDIF = 1.33-1.46*ATMTR
    elif (ATMTR <= 0.35) and (ATMTR > 0.07):
        FRDIF = 1.-2.3*(ATMTR-0.07)**2
    else:  # ATMTR <= 0.07
        FRDIF = 1.

    DIFPP = FRDIF*ATMTR*0.5*SC
    
    # Pack return values in namedtuple, add to cache and return
    astro_nt = namedtuple("AstroResults","DAYL, DAYLP, SINLD, COSLD, DIFPP, "
                                         "ATMTR, DSINBE, ANGOT")
    retvalue = astro_nt(DAYL, DAYLP, SINLD, COSLD, DIFPP, ATMTR, DSINBE, ANGOT)
    _cache[(IDAY, LAT, AVRAD)] = retvalue

    return retvalue

#-------------------------------------------------------------------------------
class Afgen(object):
    """Emulates the AFGEN function in WOFOST.
    
    :param tbl_xy: List or array of XY value pairs describing the function
        the X values should be mononically increasing.
    :param unit: The interpolated values is returned with given
        `unit <http://pypi.python.org/pypi/Unum/4.1.0>`_ assigned,
        defaults to None if Unum is not used.
    
    Returns the interpolated value provided with the 
    absicca value at which the interpolation should take place.
    
    example::
    
        >>> tbl_xy = [0,0,1,1,5,10]
        >>> f =  Afgen(tbl_xy)
        >>> f(0.5)
        0.5
        >>> f(1.5)
        2.125
        >>> f(5)
        10.0
        >>> f(6)
        10.0
        >>> f(-1)
        0.0
    """
    
    def _check_x_ascending(self, tbl_xy):
        """Checks that the x values are strictly ascending.
        
        Also truncates any trailing (0.,0.) pairs as a results of data coming
        from a CGMS database.
        """
        x_list = tbl_xy[0::2]
        y_list = tbl_xy[1::2]
        n = len(x_list)
        
        # Check if x range is ascending continuously
        rng = list(range(1, n))
        x_asc = [True if (x_list[i] > x_list[i-1]) else False for i in rng]
        
        # Check for breaks in the series where the ascending sequence stops.
        # Only 0 or 1 breaks are allowed. Use the XOR operator '^' here
        n = len(x_asc)
        sum_break = sum([1 if (x0 ^ x1) else 0 for x0,x1 in zip(x_asc, x_asc[1:])])
        if sum_break == 0:
            x = x_list
            y = y_list
        elif sum_break == 1:
            x = [x_list[0]]
            y = [y_list[0]]
            for i,p in zip(rng, x_asc):
                if p is True:
                    x.append(x_list[i])
                    y.append(y_list[i])
        else:
            msg = ("X values for AFGEN input list not strictly ascending: %s"
                   % x_list)
            raise ValueError(msg)
        
        return x, y            

    def __init__(self, tbl_xy, unit=None):
        
        self.unit = unit

        x_list, y_list = self._check_x_ascending(tbl_xy)
        x_list = self.x_list = list(map(float, x_list))
        y_list = self.y_list = list(map(float, y_list))
        intervals = list(zip(x_list, x_list[1:], y_list, y_list[1:]))
        self.slopes = [(y2 - y1)/(x2 - x1) for x1, x2, y1, y2 in intervals]

    def __call__(self, x):

        if x <= self.x_list[0]:
            return self.y_list[0]
        if x >= self.x_list[-1]:
            return self.y_list[-1]

        i = bisect_left(self.x_list, x) - 1
        v = self.y_list[i] + self.slopes[i] * (x - self.x_list[i])

        # if a unum unit is defined, multiply with a unit
        if self.unit is not None:
            v *= self.unit

        return v

#-------------------------------------------------------------------------------
class Afgen2(object):
    """Emulates the AFGEN function in TTUTIL with Numpy.interp
    
    :param tbl_xy: List or array of XY value pairs describing the function
        the X values should be mononically increasing.
    :param unit: The interpolated values is returned with given
        `unit <http://pypi.python.org/pypi/Unum/4.1.0>`_ assigned.
    
    Returns the interpolated value provided with the 
    absicca value at which the interpolation should take place.
    
    example::
    
        >>> tbl_xy = [0,0,1,1,5,10]
        >>> f =  Afgen(tbl_xy)
        >>> f(0.5)
        0.5
        >>> f(1.5)
        2.125
        >>> f(5)
        10.0
        >>> f(6)
        10.0
        >>> f(-1)
        0.0
    """
    
    def __init__(self, tbl_xy, unit=None):
    
        x = tbl_xy[0::2]
        y = tbl_xy[1::2]
        
        # Determine if there are empty pairs in tbl_xy by searching
        # for the point where x stops monotonically increasing
        xpos = np.arange(1, len(x))
        ibreak = False
        for i in xpos:
            if x[i] <= x[i-1]:
                ibreak = True
                break
        
        if ibreak is True:
            x = x[0:i]
            y = y[0:i]
        
        self.x = x
        self.y = y
        self.unit = unit
    
    def __call__(self, p):
        
        v = np.interp(p, self.x, self.y)
        
        # if a unum unit is defined, multiply with a unit
        if self.unit is not None:
            v *= self.unit
        return float(v)
    
    def __str__(self):
        msg = "AFGEN interpolation over (X,Y) pairs:\n"
        for (x,y) in zip(self.x, self.y):
            msg += ("(%f,%f)\n " % (x,y))
        msg += "\n"
        if self.unit is not None:
             msg += "Return value as unit: %s" % self.unit
        return msg

#-------------------------------------------------------------------------------
class Chainmap(MutableMapping):
    """Combine multiple mappings for sequential lookup.

    For example, to emulate Python's normal lookup sequence:

        import __builtin__
        pylookup = Chainmap(locals(), globals(), vars(__builtin__))
    
    recipe taken from http://code.activestate.com/recipes/305268/
    Available natively in the connections module in python 3.3
    """

    def __init__(self, *maps):
        self._maps = maps

    def __getitem__(self, key):
        for mapping in self._maps:
            try:
                return mapping[key]
            except KeyError:
                pass
        raise KeyError(key)
        
    def keys(self):
        k = []
        for mapping in self._maps:
            return k.extend(mapping.keys)
        return k


def merge_dict(d1, d2, overwrite=False):
    """Merge contents of d1 and d2 and return the merged dictionary
    
    Note:
    
    * The dictionaries d1 and d2 are unaltered.
    * If `overwrite=False` (default), a `RuntimeError` will be raised when
      duplicate keys exist, else any existing keys in d1 are silently
      overwritten by d2.  
    """
    # Note: May  partially be replaced by a ChainMap as of python 3.3 
    if overwrite is False:
        sd1 = set(d1.keys())
        sd2 = set(d2.keys())
        intersect = sd1.intersection(sd2)
        if len(intersect) > 0:
            msg = "Dictionaries to merge have overlapping keys: %s"
            raise RuntimeError(msg % intersect)

    td = copy.deepcopy(d1)
    td.update(d2)
    return td


class ConfigurationLoader(object):
    """Class for loading the model configuration from a PCSE configuration files

        :param config: string given file name containing model configuration
        """
    _required_attr = ("CROP", "SOIL", "AGROMANAGEMENT", "OUTPUT_VARS", "OUTPUT_INTERVAL",
                      "OUTPUT_INTERVAL_DAYS", "SUMMARY_OUTPUT_VARS")
    defined_attr = []
    model_config_file = None
    description = None

    def __init__(self, config):

        if not isinstance(config, str):
            msg = ("Keyword 'config' should provide the name of the file " +
                   "storing the configuration of the model PCSE should run.")
            raise exc.PCSEError(msg)

        # check if model configuration file is an absolute or relative path. If
        # not assume that it is located in the 'conf/' folder in the PCSE
        # distribution
        if os.path.isabs(config):
            mconf = config
        elif config.startswith("."):
            mconf = os.path.normpath(config)
        else:
            pcse_dir = os.path.dirname(__file__)
            mconf = os.path.join(pcse_dir, "conf", config)
        model_config_file = os.path.abspath(mconf)

        # check that configuration file exists
        if not os.path.exists(model_config_file):
            msg = "PCSE model configuration file does not exist: %s" % model_config_file
            raise exc.PCSEError(msg)
        # store for later use
        self.model_config_file = model_config_file

        # Load file using execfile
        try:
            loc = {}
            bytecode = compile(open(model_config_file).read(), model_config_file, 'exec')
            exec(bytecode, {}, loc)
        except Exception as e:
            msg = "Failed to load configuration from file '%s' due to: %s"
            msg = msg % (model_config_file, e)
            raise exc.PCSEError(msg)

        # Add the descriptive header for later use
        if "__doc__" in loc:
            desc = loc.pop("__doc__")
            if len(desc) > 0:
                self.description = desc
                if self.description[-1] != "\n":
                    self.description += "\n"

        # Loop through the attributes in the configuration file
        for key, value in list(loc.items()):
            if key.isupper():
                self.defined_attr.append(key)
                setattr(self, key, value)

        # Check for any missing compulsary attributes
        req = set(self._required_attr)
        diff = req.difference(set(self.defined_attr))
        if diff:
            msg = "One or more compulsary configuration items missing: %s" % list(diff)
            raise exc.PCSEError(msg)

    def __str__(self):
        msg = "PCSE ConfigurationLoader from file:\n"
        msg += "  %s\n\n" % self.model_config_file
        if self.description is not None:
            msg += ("%s Header of configuration file %s\n"% ("-"*20, "-"*20))
            msg += self.description
            if msg[-1] != "\n":
                msg += "\n"
            msg += ("%s Contents of configuration file %s\n"% ("-"*19, "-"*19))
        for k in self.defined_attr:
             r = "%s: %s" % (k, getattr(self, k))
             msg += (textwrap.fill(r, subsequent_indent="  ") + "\n")
        return msg


def is_a_month(day):
    """Returns True if the date is on the last day of a month."""

    if day.month==12:
        if day == datetime.date(day.year, day.month, 31):
            return True
    else:
        if (day == datetime.date(day.year, day.month+1, 1) - \
                   datetime.timedelta(days=1)):
            return True

    return False

def is_a_week(day, weekday=0):
    """Default weekday is Monday. Monday is 0 and Sunday is 6"""
    if day.weekday() == weekday:
        return True
    else:
        return False

def is_a_dekad(day):
    """Returns True if the date is on a dekad boundary, i.e. the 10th,
    the 20th or the last day of each month"""

    if day.month == 12:
        if day == datetime.date(day.year, day.month, 10):
            return True
        elif day == datetime.date(day.year, day.month, 20):
            return True
        elif day == datetime.date(day.year, day.month, 31):
            return True
    else:
        if day == datetime.date(day.year, day.month, 10):
            return True
        elif day == datetime.date(day.year, day.month, 20):
            return True
        elif (day == datetime.date(day.year, day.month+1, 1) -
                     datetime.timedelta(days=1)):
            return True

    return False

def load_SQLite_dump_file(dump_file_name, SQLite_db_name):
    """Build an SQLite database <SQLite_db_name> from dump file <dump_file_name>.
    """

    with open(dump_file_name) as fp:
        sql_dump = fp.readlines()
    str_sql_dump = "".join(sql_dump)
    con = sqlite3.connect(SQLite_db_name)
    con.executescript(str_sql_dump)
    con.close()

def safe_float(x):
    """Returns the value of x converted to float, if fails return None.
    """
    try:
        return float(x)
    except (ValueError, TypeError):
        return None

def check_date(indate):
        """Check representations of date and try to force into a datetime.date

        The following formats are supported:

        1. a date object
        2. a datetime object
        3. a string of the format YYYYMMDD
        4. a string of the format YYYYDDD

        Formats 2-4 are all converted into a date object internally.
        """

        import datetime as dt
        if isinstance(indate, dt.datetime):
            return indate.date()
        elif isinstance(indate, dt.date):
            return indate
        elif isinstance(indate, str):
            skey = indate.strip()
            l = len(skey)
            if l==8:
                # assume YYYYMMDD
                dkey = dt.datetime.strptime(skey,"%Y%m%d")
                return dkey.date()
            elif l==7:
                # assume YYYYDDD
                dkey = dt.datetime.strptime(skey,"%Y%j")
                return dkey.date()
            else:
                msg = "Input value not recognized as date: %s"
                raise KeyError(msg % indate)
        else:
            msg = "Input value not recognized as date: %s"
            raise KeyError(msg % indate)


class DummySoilDataProvider(dict):
    """This class is to provide some dummy soil parameters which are needed for potential production simulation
    """
    _defaults = {"SMFCF":0.3,
                 "SM0":0.4,
                 "SMW":0.1,
                 "RDMSOL":120,
                 "CRAIRC":0.06,
                 "K0":10.,
                 "SOPE":10.,
                 "KSUB":10.}

    def __init__(self):
        self.update(self._defaults)