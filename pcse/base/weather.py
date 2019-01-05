# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Base classes for creating PCSE simulation units.

In general these classes are not to be used directly, but are to be subclassed
when creating PCSE simulation units.
"""
import logging
import pickle

from .. import exceptions as exc
from ..settings import settings


class SlotPickleMixin(object):
    """This mixin makes it possible to pickle/unpickle objects with __slots__ defined.

    In many programs, one or a few classes have a very large number of instances.
    Adding __slots__ to these classes can dramatically reduce the memory footprint
    and improve execution speed by eliminating the instance dictionary. Unfortunately,
    the resulting objects cannot be pickled. This mixin makes such classes pickleable
    again and even maintains compatibility with pickle files created before adding
    __slots__.

    Recipe taken from:
    http://code.activestate.com/recipes/578433-mixin-for-pickling-objects-with-__slots__/
    """

    def __getstate__(self):
        return dict(
            (slot, getattr(self, slot))
            for slot in self.__slots__
            if hasattr(self, slot)
        )

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)


class WeatherDataContainer(SlotPickleMixin):
    """Class for storing weather data elements.

    Weather data elements are provided through keywords that are also the
    attribute names under which the variables can accessed in the
    WeatherDataContainer. So the keyword TMAX=15 sets an attribute
    TMAX with value 15.

    The following keywords are compulsory:

    :keyword LAT: Latitude of location (decimal degree)
    :keyword LON: Longitude of location (decimal degree)
    :keyword ELEV: Elevation of location (meters)
    :keyword DAY: the day of observation (python datetime.date)
    :keyword IRRAD: Incoming global radiaiton (J/m2/day)
    :keyword TMIN: Daily minimum temperature (Celsius)
    :keyword TMAX: Daily maximum temperature (Celsius)
    :keyword VAP: Daily mean vapour pressure (hPa)
    :keyword RAIN: Daily total rainfall (cm/day)
    :keyword WIND: Daily mean wind speed at 2m height (m/sec)
    :keyword E0: Daily evaporation rate from open water (cm/day)
    :keyword ES0: Daily evaporation rate from bare soil (cm/day)
    :keyword ET0: Daily evapotranspiration rate from reference crop (cm/day)

    There are two optional keywords arguments:

    :keyword TEMP: Daily mean temperature (Celsius), will otherwise be
                   derived from (TMAX+TMIN)/2.
    :keyword SNOWDEPTH: Depth of snow cover (cm)
    """
    sitevar = ["LAT", "LON", "ELEV"]
    required = ["IRRAD", "TMIN", "TMAX", "VAP", "RAIN", "E0", "ES0", "ET0", "WIND"]
    optional = ["SNOWDEPTH", "TEMP", "TMINRA"]
    # In the future __slots__ can be extended or attribute setting can be allowed
    # by add '__dict__' to __slots__.
    __slots__ = sitevar + required + optional + ["DAY"]

    units = {"IRRAD": "J/m2/day", "TMIN": "Celsius", "TMAX": "Celsius", "VAP": "hPa",
             "RAIN": "cm/day", "E0": "cm/day", "ES0": "cm/day", "ET0": "cm/day",
             "LAT": "Degrees", "LON": "Degrees", "ELEV": "m", "SNOWDEPTH": "cm",
             "TEMP": "Celsius", "TMINRA": "Celsius", "WIND": "m/sec"}

    # ranges for meteorological variables
    ranges = {"LAT": (-90., 90.),
              "LON": (-180., 180.),
              "ELEV": (-300, 6000),
              "IRRAD": (0., 40e6),
              "TMIN": (-50., 60.),
              "TMAX": (-50., 60.),
              "VAP": (0.06, 199.3),  # hPa, computed as sat. vapour pressure at -50, 60 Celsius
              "RAIN": (0, 25),
              "E0": (0., 2.5),
              "ES0": (0., 2.5),
              "ET0": (0., 2.5),
              "WIND": (0., 100.),
              "SNOWDEPTH": (0., 250.),
              "TEMP": (-50., 60.),
              "TMINRA": (-50., 60.)}

    def __init__(self, *args, **kwargs):

        # only keyword parameters should be used for weather data container
        if len(args) > 0:
            msg = ("WeatherDataContainer should be initialized by providing weather " +
                   "variables through keywords only. Got '%s' instead.")
            raise exc.PCSEError(msg % args)

        # First assign site variables
        for varname in self.sitevar:
            try:
                setattr(self, varname, float(kwargs.pop(varname)))
            except (KeyError, ValueError) as e:
                msg = "Site parameter '%s' missing or invalid when building WeatherDataContainer: %s"
                raise exc.PCSEError(msg, varname, e)

        # check if we have a DAY element
        if "DAY" not in kwargs:
            msg = "Date of observations 'DAY' not provided when building WeatherDataContainer."
            raise exc.PCSEError(msg)
        self.DAY = kwargs.pop("DAY")

        # Loop over required arguments to see if all required variables are there
        for varname in self.required:
            value = kwargs.pop(varname, None)
            try:
                setattr(self, varname, float(value))
            except (KeyError, ValueError, TypeError) as e:
                msg = "%s: Weather attribute '%s' missing or invalid numerical value: %s"
                logging.warning(msg, self.DAY, varname, value)

        # Loop over optional arguments
        for varname in self.optional:
            value = kwargs.pop(varname, None)
            if value is None:
                continue
            else:
                try:
                    setattr(self, varname, float(value))
                except (KeyError, ValueError, TypeError) as e:
                    msg = "%s: Weather attribute '%s' missing or invalid numerical value: %s"
                    logging.warning(msg, self.DAY, varname, value)

        # Check for remaining unknown arguments
        if len(kwargs) > 0:
            msg = "WeatherDataContainer: unknown keywords '%s' are ignored!"
            logging.warning(msg, kwargs.keys())

    def __setattr__(self, key, value):
        # Override to allow range checking on known meteo variables.

        # Skip range checking if disabled by user
        if settings.METEO_RANGE_CHECKS:
            if key in self.ranges:
                vmin, vmax = self.ranges[key]
                if not vmin <= value <= vmax:
                    msg = "Value (%s) for meteo variable '%s' outside allowed range (%s, %s)." % (
                    value, key, vmin, vmax)
                    raise exc.PCSEError(msg)
        SlotPickleMixin.__setattr__(self, key, value)

    def __str__(self):
        msg = "Weather data for %s (DAY)\n" % self.DAY
        for v in self.required:
            value = getattr(self, v, None)
            if value is None:
                msg += "%5s: element missing!\n"
            else:
                unit = self.units[v]
                msg += "%5s: %12.2f %9s\n" % (v, value, unit)
        for v in self.optional:
            value = getattr(self, v, None)
            if value is None:
                continue
            else:
                unit = self.units[v]
                msg += "%5s: %12.2f %9s\n" % (v, value, unit)
        msg += ("Latitude  (LAT): %8.2f degr.\n" % self.LAT)
        msg += ("Longitude (LON): %8.2f degr.\n" % self.LON)
        msg += ("Elevation (ELEV): %6.1f m.\n" % self.ELEV)
        return msg

    def add_variable(self, varname, value, unit):
        """Adds an attribute <varname> with <value> and given <unit>

        :param varname: Name of variable to be set as attribute name (string)
        :param value: value of variable (attribute) to be added.
        :param unit: string representation of the unit of the variable. Is
            only use for printing the contents of the WeatherDataContainer.
        """
        if varname not in self.units:
            self.units[varname] = unit
        setattr(self, varname, value)


class WeatherDataProvider(object):
    """Base class for all weather data providers.

    Support for weather ensembles in a WeatherDataProvider has to be indicated
    by setting the class variable `supports_ensembles = True`

    Example::

        class MyWeatherDataProviderWithEnsembles(WeatherDataProvider):
            supports_ensembles = True

            def __init__(self):
                WeatherDataProvider.__init__(self)

                # remaining initialization stuff goes here.
    """
    supports_ensembles = False

    # Descriptive items for a WeatherDataProvider
    longitude = None
    latitude = None
    elevation = None
    description = []
    _first_date = None
    _last_date = None
    angstA = None
    angstB = None
    # model used for reference ET
    ETmodel = "PM"

    def __init__(self):
        self.store = {}

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)

    def _dump(self, cache_fname):
        """Dumps the contents into cache_fname using pickle.

        Dumps the values of self.store, longitude, latitude, elevation and description
        """
        with open(cache_fname, "wb") as fp:
            dmp = (self.store, self.elevation, self.longitude, self.latitude, self.description, self.ETmodel)
            pickle.dump(dmp, fp, pickle.HIGHEST_PROTOCOL)

    def _load(self, cache_fname):
        """Loads the contents from cache_fname using pickle.

        Loads the values of self.store, longitude, latitude, elevation and description
        from cache_fname and also sets the self.first_date, self.last_date
        """

        with open(cache_fname, "rb") as fp:
            (store, self.elevation, self.longitude, self.latitude, self.description, ETModel) = pickle.load(fp)

        # Check if the reference ET from the cache file is calculated with the same model as
        # specified by self.ETmodel
        if ETModel != self.ETmodel:
            msg = "Mismatch in reference ET from cache file."
            raise exc.PCSEError(msg)

        self.store.update(store)

    def export(self):
        """Exports the contents of the WeatherDataProvider as a list of dictionaries.

        The results from export can be directly converted to a Pandas dataframe
        which is convenient for plotting or analyses.
        """
        weather_data = []
        if self.supports_ensembles:
            # We have to include the member_id in each dict with weather data
            pass
        else:
            days = sorted([r[0] for r in self.store.keys()])
            for day in days:
                wdc = self(day)
                r = {key: getattr(wdc, key) for key in wdc.__slots__ if hasattr(wdc, key)}
                weather_data.append(r)
        return weather_data

    @property
    def first_date(self):
        try:
            self._first_date = min(self.store)[0]
        except ValueError:
            pass
        return self._first_date

    @property
    def last_date(self):
        try:
            self._last_date = max(self.store)[0]
        except ValueError:
            pass
        return self._last_date

    @property
    def missing(self):
        missing = (self.last_date - self.first_date).days - len(self.store) + 1
        return missing

    def check_keydate(self, key):
        """Check representations of date for storage/retrieval of weather data.

        The following formats are supported:

        1. a date object
        2. a datetime object
        3. a string of the format YYYYMMDD
        4. a string of the format YYYYDDD

        Formats 2-4 are all converted into a date object internally.
        """

        import datetime as dt
        if isinstance(key, dt.datetime):
            return key.date()
        elif isinstance(key, dt.date):
            return key
        elif isinstance(key, (str, int)):
            skey = str(key).strip()
            l = len(skey)
            if l == 8:
                # assume YYYYMMDD
                dkey = dt.datetime.strptime(skey, "%Y%m%d")
                return dkey.date()
            elif l == 7:
                # assume YYYYDDD
                dkey = dt.datetime.strptime(skey, "%Y%j")
                return dkey.date()
            else:
                msg = "Key for WeatherDataProvider not recognized as date: %s"
                raise KeyError(msg % key)
        else:
            msg = "Key for WeatherDataProvider not recognized as date: %s"
            raise KeyError(msg % key)

    def _store_WeatherDataContainer(self, wdc, keydate, member_id=0):
        """Stores the WDC under given keydate and member_id.
        """

        if member_id != 0 and self.supports_ensembles is False:
            msg = "Storing ensemble weather is not supported."
            raise exc.WeatherDataProviderError(msg)

        kd = self.check_keydate(keydate)
        if not (isinstance(member_id, int) and member_id >= 0):
            msg = "Member id should be a positive integer, found %s" % member_id
            raise exc.WeatherDataProviderError(msg)

        self.store[(kd, member_id)] = wdc

    def __call__(self, day, member_id=0):

        if self.supports_ensembles is False and member_id != 0:
            msg = "Retrieving ensemble weather is not supported by %s" % self.__class__.__name__
            raise exc.WeatherDataProviderError(msg)

        keydate = self.check_keydate(day)
        if self.supports_ensembles is False:
            msg = "Retrieving weather data for day %s" % keydate
            self.logger.debug(msg)
            try:
                return self.store[(keydate, 0)]
            except KeyError as e:
                msg = "No weather data for %s." % keydate
                raise exc.WeatherDataProviderError(msg)
        else:
            msg = "Retrieving ensemble weather data for day %s member %i" % \
                  (keydate, member_id)
            self.logger.debug(msg)
            try:
                return self.store[(keydate, member_id)]
            except KeyError:
                msg = "No weather data for (%s, %i)." % (keydate, member_id)
                raise exc.WeatherDataProviderError(msg)

    def __str__(self):

        msg = "Weather data provided by: %s\n" % self.__class__.__name__
        msg += "--------Description---------\n"
        if isinstance(self.description, str):
            msg += ("%s\n" % self.description)
        else:
            for l in self.description:
                msg += ("%s\n" % str(l))
        msg += "----Site characteristics----\n"
        msg += "Elevation: %6.1f\n" % self.elevation
        msg += "Latitude:  %6.3f\n" % self.latitude
        msg += "Longitude: %6.3f\n" % self.longitude
        msg += "Data available for %s - %s\n" % (self.first_date, self.last_date)
        msg += "Number of missing days: %i\n" % self.missing
        return msg

