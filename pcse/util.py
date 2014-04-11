"""Miscelaneous utilities for PyWOFOST
"""
import os, sys
import datetime
import copy
from math import log10, cos, sin, asin, sqrt, exp
from collections import namedtuple
from bisect import bisect_left
import UserDict
import textwrap

import numpy as np
from sqlalchemy import Column

from . import exceptions as exc

def penman(day, LAT, ELEV, ANGSTA, ANGSTB, TMIN, TMAX, AVRAD, VAP, WIND2):
    """Calculates E0, ES0, ET0 based on the Penman model.
    
     This routine calculates the potential evapo(transpi)ration rates from
     a free water surface (E0), a bare soil surface (ES0), and a crop canopy
     (ET0) in mm/d. For these calculations the analysis by Penman is followed
     (Frere and Popov, 1979;Penman, 1948, 1956, and 1963).
     Subroutines and functions called: ASTRO, LIMIT.

    Input variables::
    
        day     I4  Python date                                    -      
        LAT     R4  Latitude of the site                        degrees   
        ELEV    R4  Elevation above sea level                      m      
        ANGSTA  R4  Empirical constant in Angstrom formula         -      
        ANGSTB  R4  Empirical constant in Angstrom formula         -      
        TMIN    R4  Minimum temperature                            C      
        TMAX    R4  Maximum temperature                            C      
        AVRAD   R4  Daily shortwave radiation                   J m-2 d-1 
        VAP     R4  24 hour average vapour pressure               hPa     
        WIND2   R4  24 hour average windspeed at 2 meter          m/s     
        
    Output is a tuple (E0,ES0,ET0)::
    
        E0      R4  Penman potential evaporation from a free water surface [mm/d]
        ES0     R4  Penman potential evaporation from a moist bare soil surface [mm/d]
        ET0     R4  Penman potential transpiration from a crop canopy [mm/d]
    """
    # psychrometric instrument constant (mbar/Celsius-1)
    # albedo for water surface, soil surface and canopy
    # latent heat of evaporation of water (J/kg=J/mm)
    # Stefan Boltzmann constant (J/m2/d/K4)
    PSYCON=0.67; REFCFW=0.05; REFCFS=0.15; REFCFC=0.25
    LHVAP=2.45E6; STBC=4.9E-3

    # preparatory calculations
    # mean daily temperature and temperature difference (Celsius)
    # coefficient Bu in wind function, dependent on temperature
    # difference
    TMPA  = (TMIN+TMAX)/2.
    TDIF  = TMAX - TMIN
    BU    = 0.54 + 0.35 * limit(0.,1.,(TDIF-12.)/4.)

    # barometric pressure (mbar)
    # psychrometric constant (mbar/Celsius)
    PBAR  = 1013.*exp(-0.034*ELEV/(TMPA+273.))
    GAMMA = PSYCON*PBAR/1013.


    # saturated vapour pressure according to equation of Goudriaan
    # (1977) derivative of SVAP with respect to temperature, i.e. 
    # slope of the SVAP-temperature curve (mbar/Celsius);
    # measured vapour pressure not to exceed saturated vapour pressure

    SVAP  = 6.10588 * exp(17.32491*TMPA/(TMPA+238.102))
    DELTA = 238.102*17.32491*SVAP/(TMPA+238.102)**2
    VAP   = min(VAP,SVAP)

    # the expression n/N (RELSSD) from the Penman formula is estimated
    # from the Angstrom formula: RI=RA(A+B.n/N) -> n/N=(RI/RA-A)/B,
    # where RI/RA is the atmospheric transmission obtained by a CALL
    # to ASTRO:

    r = astro(day, LAT, AVRAD)
    RELSSD = limit(0.,1.,(r.ATMTR-abs(ANGSTA))/abs(ANGSTB))

    # Terms in Penman formula, for water, soil and canopy

    # net outgoing long-wave radiation (J/m2/d) acc. to Brunt (1932)
    RB  = STBC*(TMPA+273.)**4*(0.56-0.079*sqrt(VAP))*(0.1+0.9*RELSSD)

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
    
    return (E0,ES0,ET0)

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
    if (A < MIN_A or A > MAX_A):
        msg = "invalid Angstrom A value!"
        raise exc.PCSEError(msg)
    if (B < MIN_B or B > MAX_B):
        msg = "invalid Angstrom B value!"
        raise exc.PCSEError(msg)
    if (SUM_AB < MIN_SUM_AB or SUM_AB > MAX_SUM_AB):
        msg = "invalid sum of Angstrom A & B values!"
        raise exc.PCSEError(msg)

    return [A,B]

#-------------------------------------------------------------------------------
def wind10to2(wind10):
    """Converts windspeed at 10m to windspeed at 2m using log. wind profile
    """
    wind2 = wind10 * (log10(2./0.033) / log10(10/0.033))
    return wind2

#-------------------------------------------------------------------------------
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
        raise ValueError, 'tdew=%g is not in range -95 to +60 deg C' % tdew

    tmp = (17.27 * tdew) / (tdew + 237.3)
    ea = 0.6108 * exp(tmp)
    return ea
#-------------------------------------------------------------------------------
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
                       
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
def is_data_column(col):
    """Returns True if given column is not part of the run descriptors
    
    Parameters:
    * col - An SQLAlchemy column object
    """
    if not isinstance(col, Column):
        msg = "Provided variable is not an SQLAlchemy column object"
        raise RuntimeError(msg)

    run_descriptors = ["grid_no","crop_no","year","day","simulation_mode", 
                       "member_id"]
    if col.name not in run_descriptors:
        return True
    else:
        return False

#-------------------------------------------------------------------------------
def limit(min, max, v):
    """limits the range of v between min and max
    """

    if min > max:
        raise RuntimeError("Min value larger than max")
    
    if v < min:       # V below range: return min
        return min
    elif v < max:     # v within range: return v
        return v
    else:             # v above range: return max
        return max

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
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
    if (abs(AOB) <= 1.0):
        DAYL  = 12.0*(1.+2.*ASIN(AOB)/PI)
        # integrals of sine of solar height
        DSINB  = 3600.*(DAYL*SINLD+24.*COSLD*SQRT(1.-AOB**2)/PI)
        DSINBE = 3600.*(DAYL*(SINLD+0.4*(SINLD**2+COSLD**2*0.5))+
                 12.*COSLD*(2.+3.*0.4*SINLD)*SQRT(1.-AOB**2)/PI)
    else:
        if (AOB >  1.0): DAYL = 24.0
        if (AOB < -1.0): DAYL =  0.0
        # integrals of sine of solar height	
        DSINB  = 3600.*(DAYL*SINLD)
        DSINBE = 3600.*(DAYL*(SINLD+0.4*(SINLD**2+COSLD**2*0.5)))

    # Calculate solution for base=-4 (ANGLE) degrees
    AOB_CORR = (-SIN(ANGLE*RAD)+SINLD)/COSLD
    if (ABS(AOB_CORR) <= 1.0):
        DAYLP = 12.0*(1.+2.*ASIN(AOB_CORR)/PI)
    elif (AOB_CORR >  1.0):
        DAYLP = 24.0
    elif (AOB_CORR < -1.0):
        DAYLP =  0.0

    # extraterrestrial radiation and atmospheric transmission
    ANGOT  = SC*DSINB
    # Check for DAYL=0 as in that case the angot radiation is 0 as well
    if (DAYL > 0.0):
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
    elif (ATMTR <= 0.07):
        FRDIF = 1.

    DIFPP = FRDIF*ATMTR*0.5*SC
    
    # Pack return values in namedtuple, add to cache and return
    astro_nt = namedtuple("AstroResults","DAYL,DAYLP,SINLD,COSLD,DIFPP,ATMTR,DSINBE,ANGOT")
    retvalue = astro_nt(DAYL,DAYLP,SINLD,COSLD,DIFPP,ATMTR,DSINBE,ANGOT)
    _cache[(IDAY, LAT, AVRAD)] = retvalue

    #print IDAY, LAT, AVRAD, DAYL, DAYLP, SINLD, COSLD, DIFPP, ATMTR, DSINBE

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
        rng = range(1, n)
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
        x_list = self.x_list = map(float, x_list)
        y_list = self.y_list = map(float, y_list)
        intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
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
        xpos = np.arange(len(x)-1) + 1
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
class Chainmap(UserDict.DictMixin):
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

#-------------------------------------------------------------------------------
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
            execfile(model_config_file, {}, loc)
        except Exception, e:
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
        for key, value in loc.items():
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
        if (day == datetime.date(day.year, day.month, 31)):
            return True
    else:
        if (day == datetime.date(day.year, day.month+1, 1) - \
                   datetime.timedelta(days=1)):
            return True

    return False


def is_a_dekad(day):
    """Returns True if the date is on a dekad boundary, i.e. the 10th,
    the 20th or the last day of each month"""

    if day.month==12:
        if (day == datetime.date(day.year, day.month, 10)):
            return True
        elif (day == datetime.date(day.year, day.month, 20)):
            return True
        elif (day == datetime.date(day.year, day.month, 31)):
            return True
    else:
        if (day == datetime.date(day.year, day.month, 10)):
            return True
        elif (day == datetime.date(day.year, day.month, 20)):
            return True
        elif (day == datetime.date(day.year, day.month+1, 1) - \
                   datetime.timedelta(days=1)):
            return True

    return False
