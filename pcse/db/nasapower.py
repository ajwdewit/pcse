# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import sys, os
import urllib
from datetime import date, timedelta
import numpy as np

from ..base_classes import WeatherDataProvider, WeatherDataContainer
from ..util import wind10to2, ea_from_tdew, reference_ET, check_angstromAB
from ..exceptions import PCSEError
from ..settings import settings

# Define some lambdas to take care of unit conversions.
MJ_to_J = lambda x: x * 1e6
no_conv = lambda x: x
mm_to_cm = lambda x: x / 10.
tdew_to_hpa = lambda x: ea_from_tdew(x) * 10.


class NASAPowerWeatherDataProvider(WeatherDataProvider):
    """WeatherDataProvider for using the NASA POWER database with PCSE

    :param latitude: latitude to request weather data for
    :param longitude: longitude to request weather data for
    :keyword force_update: Set to True to force to request fresh data
        from POWER website.
    :keyword ETmodel: "PM"|"P" for selecting penman-monteith or Penman
        method for reference evapotranspiration. Defaults to "PM".

    The NASA POWER database is a global database of daily weather data
    specifically designed for agrometeorological applications. The spatial
    resolution of the database is 1x1 degrees. It is derived from weather
    station observations in combination with satellite data for parameters
    like radiation.

    The weather data is updated with a delay of about 3 months which makes
    the database unsuitable for real-time monitoring, nevertheless the
    POWER database is useful for many other studies and it is a major
    improvement compared to the monthly weather data that were used with
    WOFOST in the past.

    For more information on the NASA POWER database see the documentation
    at: http://power.larc.nasa.gov/common/AgroclimatologyMethodology/Agro_Methodology_Content.html

    The `NASAPowerWeatherDataProvider` retrieves the weather from the
    th NASA POWER website and does the necessary conversions to be compatible
    with PCSE. After the data has been retrieved and stored, the contents
    are dumped to a binary cache file. If another request is made for the
    same location, the cache file is loaded instead of a full request to the
    NASA Power server.

    Cache files are used until they are older then 90 days. After 90 days
    the NASAPowerWeatherDataProvider will make a new request to obtain
    more recent data from the NASA POWER server. If this request fails
    it will fall back to the existing cache file. The update of the cache
    file can be forced by setting `force_update=True`.

    Finally, note that any latitude/longitude within a 1x1 degrees grid box
    will yield the same weather data, e.g. there is no difference between
    lat/lon 5.3/52.1 and lat/lon 5.9/52.8. Nevertheless slight differences
    in PCSE simulations may occur due to small differences in day length.

    """
    # Variable names in POWER data
    power_variables = ["toa_dwn", "swv_dwn", "lwv_dwn", "T2M", "T2MN",
                       "T2MX", "RH2M", "DFP2M", "RAIN", "WS10M"]

    # Mapping PCSE name to power name, conversion factor and unit of weather variables
    pcse_variables = [("IRRAD", "swv_dwn", MJ_to_J, "J/m2/day"),
                      ("TMIN", "T2MN", no_conv, "Celsius"),
                      ("TMAX", "T2MX", no_conv, "Celsius"),
                      ("TEMP", "T2M", no_conv, "Celsius"),
                      ("VAP", "DFP2M", tdew_to_hpa, "hPa"),
                      ("WIND", "WS10M", wind10to2, "m/sec"),
                      ("RAIN", "RAIN", mm_to_cm, "cm/day")]
    # other constants
    HTTP_OK = 200
    angstA = 0.25
    angstB = 0.5

    def __init__(self, latitude, longitude, force_update=False, ETmodel="PM"):

        WeatherDataProvider.__init__(self)

        if latitude < -90 or latitude > 90:
            msg = "Latitude should be between -90 and 90 degrees."
            raise ValueError(msg)
        if longitude < -180 or longitude > 180:
            msg = "Longitude should be between -180 and 180 degrees."
            raise ValueError(msg)

        self.latitude = float(latitude)
        self.longitude = float(longitude)
        self.ETmodel = ETmodel
        msg = "Retrieving weather data from NASA Power for lat/lon: (%f, %f)."
        self.logger.info(msg % (self.latitude, self.longitude))

        # Check for existence of a cache file
        cache_file = self._find_cache_file(self.latitude, self.longitude)
        if cache_file is None or force_update is True:
            msg = "No cache file or forced update, getting data from NASA Power."
            self.logger.debug(msg)
            # No cache file, we really have to get the data from the NASA server
            self._get_and_process_NASAPower(self.latitude, self.longitude)
            return

        # get age of cache file, if age < 90 days then try to load it. If loading fails retrieve data
        # from the NASA server .
        r = os.stat(cache_file)
        cache_file_date = date.fromtimestamp(r.st_mtime)
        age = (date.today() - cache_file_date).days
        if age < 90:
            msg = "Start loading weather data from cache file: %s" % cache_file
            self.logger.debug(msg)

            status = self._load_cache_file()
            if status is not True:
                msg = "Loading cache file failed, reloading data from NASA Power."
                self.logger.debug(msg)
                # Loading cache file failed!
                self._get_and_process_NASAPower(self.latitude, self.longitude)
        else:
            # Cache file is too old. Try loading new data from NASA
            try:
                msg = "Cache file older then 90 days, reloading data from NASA Power."
                self.logger.debug(msg)
                self._get_and_process_NASAPower(self.latitude, self.longitude)
            except Exception as e:
                msg = ("Reloading data from NASA failed, reverting to (outdated) " +
                       "cache file")
                self.logger.debug(msg)
                status = self._load_cache_file()
                if status is not True:
                    msg = "Outdated cache file failed loading."
                    raise PCSEError(msg)

    def _get_and_process_NASAPower(self, latitude, longitude):
        """Handles the retrieval and processing of the NASA Power data
        """
        powerdata = self._query_NASAPower_server(latitude, longitude)

        # Store the informational header then parse variables
        self.description = self._parse_header(powerdata)
        self.elevation = self._parse_elevation(powerdata)
        recs = self._process_power_records(powerdata)

        # Determine Angstrom A/B parameters
        self.AngstA, self.AngstB = self._estimate_AngstAB(recs)

        # Start building the weather data containers
        self._make_WeatherDataContainers(recs)

        # dump contents to a cache file
        cache_filename = self._get_cache_filename(latitude, longitude)
        self._dump(cache_filename)

    def _estimate_AngstAB(self, power_records):
        """Determine Angstrom A/B parameters from Top-of-Atmosphere (toa_dwn) and
        top-of-Canopy (swv_dwn) radiation values.

        :param power_records: rows of parsed power records (see self.power_variables)
        :return: tuple of Angstrom A/B values

        The Angstrom A/B parameters are determined by dividing swv_dwn by toa_dwn
        and taking the 0.05 percentile for Angstrom A and the 0.98 percentile for
        Angstrom A+B: toa_dwn*(A+B) approaches the upper envelope while
        toa_dwn*A approaches the lower envelope of the records of swv_dwn
        values.
        """

        msg = "Start estimation of Angstrom A/B values from POWER data."
        self.logger.debug(msg)

        # check if sufficient data is available to make a reasonable estimate:
        # As a rule of thumb we want to have at least 200 days available
        if len(power_records) < 200:
            msg = ("Less then 200 days of data available. Reverting to " +
                   "default Angstrom A/B coefficients (%f, %f)")
            self.logger.warn(msg % (self.angstA, self.angstB))
            return self.angstA, self.angstB

        # calculate relative radiation (swv_dwn/toa_dwn) and percentiles
        relative_radiation = np.array([rec["swv_dwn"]/rec["toa_dwn"] for rec in power_records])
        angstrom_a = float(np.percentile(relative_radiation, 5))
        angstrom_ab = float(np.percentile(relative_radiation, 98))
        angstrom_b = angstrom_ab - angstrom_a

        try:
            check_angstromAB(angstrom_a, angstrom_b)
        except PCSEError as e:
            msg = ("Angstrom A/B values (%f, %f) outside valid range: %s. " +
                   "Reverting to default values.")
            msg = msg % (angstrom_a, angstrom_b, e)
            self.logger.warn(msg)
            return self.angstA, self.angstB

        msg = "Angstrom A/B values estimated: (%f, %f)." % (angstrom_a, angstrom_b)
        self.logger.debug(msg)

        return angstrom_a, angstrom_b

    def _query_NASAPower_server(self, latitude, longitude):
        """Query the NASA Power server for data on given latitude/longitude
        """

        # build URL for retrieving data
        server = "power.larc.nasa.gov"
        t_url = ("http://{server}/cgi-bin/cgiwrap/solar/agro.cgi?" +
                 "email=agroclim%40larc.nasa.gov&step=1&lat={lat}&lon={lon}" +
                 "&ms=1&ds=1&ys=1984&me={month}&de={day}&ye={year}&p=toa_dwn&" +
                 "p=swv_dwn&p=lwv_dwn&p=T2M&p=T2MN&p=T2MX&p=RH2M&" +
                 "p=DFP2M&p=RAIN&p=WS10M&submit=Submit")
        d = date.today()
        url = t_url.format(server=server, lat=latitude, lon=longitude, year=d.year,
                           month=d.month, day=d.day)
        msg = "Starting retrieval from NASA Power with URL: %s" % url
        self.logger.debug(msg)
        req = urllib.urlopen(url)

        if req.getcode() != self.HTTP_OK:
            msg = ("Failed retrieving POWER data from %s. Server returned HTTP " +
                   "code: %i") % (server, req.getcode())
            raise PCSEError(msg)

        powerdata = req.readlines()
        req.close()

        msg = "Successfully retrieved data from NASA Power"
        self.logger.debug(msg)

        return powerdata

    def _find_cache_file(self, latitude, longitude):
        """Try to find a cache file for given latitude/longitude.

        Returns None if the cache file does not exist, else it returns the full path
        to the cache file.
        """
        cache_filename = self._get_cache_filename(latitude, longitude)
        if os.path.exists(cache_filename):
            return cache_filename
        else:
            return None

    def _get_cache_filename(self, latitude, longitude):
        """Constructs the filename used for cache files given latitude and longitude

        The latitude and longitude is coded into the filename by truncating on
        0.1 degree. So the cache filename for a point with lat/lon 52.56/-124.78 will be:
        NASAPowerWeatherDataProvider_LAT00525_LON-1247.cache
        """

        fname = "%s_LAT%05i_LON%05i.cache" % (self.__class__.__name__,
                                              int(latitude*10), int(longitude*10))
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, fname)
        return cache_filename

    def _write_cache_file(self):
        """Writes the meteo data from NASA Power to a cache file.
        """
        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        try:
            self._dump(cache_filename)
        except (IOError, EnvironmentError) as e:
            msg = "Failed to write cache to file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)

    def _load_cache_file(self):
        """Loads the data from the cache file. Return True if successful.
        """
        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        try:
            self._load(cache_filename)
            msg = "Cache file successfully loaded."
            self.logger.debug(msg)
            return True
        except (IOError, EnvironmentError, EOFError) as e:
            msg = "Failed to load cache from file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
            return False


    def _make_WeatherDataContainers(self, recs):
        """Create a WeatherDataContainers from recs, compute ET and store the WDC's.
        """

        for rec in recs:
            t = {"LAT": self.latitude, "LON": self.longitude, "ELEV": self.elevation, "DAY": rec["day"]}

            for pcse_name, power_name, conv, unit in self.pcse_variables:
                value = conv(rec[power_name])
                t[pcse_name] = value

            # Reference evapotranspiration in mm/day
            try:
                E0, ES0, ET0 = reference_ET(t["DAY"], t["LAT"], t["ELEV"], t["TMIN"], t["TMAX"], t["IRRAD"],
                                            t["VAP"], t["WIND"], self.angstA, self.angstB, self.ETmodel)
            except ValueError as e:
                msg = (("Failed to calculate reference ET values on %s. " % t["DAY"]) +
                       ("With input values:\n %s.\n" % str(t)) +
                       ("Due to error: %s" % e))
                raise PCSEError(msg)

            # update record with ET values value convert to cm/day
            t.update({"E0": E0/10., "ES0": ES0/10., "ET0": ET0/10.})

            # Build weather data container from dict 't'
            wdc = WeatherDataContainer(**t)

            # add wdc to dictionary for thisdate
            self._store_WeatherDataContainer(wdc, wdc.DAY)

    def _parse_elevation(self, powerdata):
        """Parse elevation out of the powerdata header"""
        for line in powerdata:
            if line.startswith(b"Elevation"):
                elev = int(line.split(b"=")[-1])
                return elev

    def _parse_header(self, powerdata):
        header = powerdata[1:7]
        header.append(powerdata[20])
        return [h.strip() for h in header]

    def _process_power_records(self, powerdata):
        """Process the meteorological records returned by NASA POWER
        """
        msg = "Start parsing of POWER records from URL retrieval."
        self.logger.debug(msg)

        # First strip off the header by searching for '-END HEADER'
        is_header = True
        while is_header:
            line = powerdata.pop(0)
            if line.startswith(b"-END HEADER"):
                is_header = False

        # Now start parsing meteo records
        recs = []
        for line in powerdata:
            rec = self._parse_raw_power_record(line)
            if rec is not None:
                recs.append(rec)

        nrecs = len(powerdata)
        nsuccess = len(recs)
        nfail = nrecs - nsuccess
        msg = "Parsed %i POWER records: %i success, %i failures" % (nrecs, nsuccess, nfail)
        self.logger.debug(msg)

        return recs

    def _parse_raw_power_record(self, line):
        """Parse each record and return as a dict, return None if parsing fails due to a missing variable
        """
        r = line.split()
        rec = {}
        year = int(r.pop(0))
        doy = int(r.pop(0))
        rec["day"] = date(year, 1, 1) + timedelta(days=(doy - 1))
        try:
            for i, name in enumerate(self.power_variables):
                rec[name] = float(r[i])
            return rec
        except ValueError:
            msg = "POWER record for day '%s' failed to parse (probably incomplete)" % rec["day"]
            self.logger.debug(msg)
            return None


