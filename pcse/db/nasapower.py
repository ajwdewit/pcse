import sys, os
import urllib
from datetime import date, timedelta
import time

from ..base_classes import WeatherDataProvider, WeatherDataContainer
from ..util import wind10to2, ea_from_tdew, penman
from ..exceptions import PCSEError
from .. import global_settings

# Define some lambdas to take care of unit conversions.
MJ_to_J = lambda x: x * 1e6
no_conv = lambda x: x
mm_to_cm = lambda x: x / 10.
tdew_to_hpa = lambda x: ea_from_tdew(x) * 10


class NASAPowerWeatherDataProvider(WeatherDataProvider):
    """
    """

    # Variable names in POWER data
    power_variables = ["toa_dwn", "swv_dwn", "lwv_dwn", "T2M", "T2MN", "T2MX", "RH2M", "DFP2M", "RAIN", "WS10M"]

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

    def __init__(self, latitude=None, longitude=None, cache_path=None,
                 force_update=False):
        """
        """
        WeatherDataProvider.__init__(self)

        # Warning for using default Angstrom A/B values
        if (self.angstA, self.angstB) == (0.25, 0.5):
            self.logger.warning("Warning: using default values for Angstrom A/B (0.25, 0.5)!")

        if latitude < -90 or latitude > 90:
            msg = "Latitude should be between -90 and 90 degrees."
            raise ValueError(msg)
        if longitude < -180 or longitude > 180:
            msg = "Longitude should be between -180 and 180 degrees."
            raise ValueError(msg)

        self.latitude = float(latitude)
        self.longitude = float(longitude)

        # Check for existence of a cache file
        cache_file = self._find_cache_file(self.latitude, self.longitude)
        if cache_file is None:
            # No cache file, we really have to get the data from the NASA server
            self._get_and_process_NASAPower(self.latitude, self.longitude)
            return

        # get age of cache file, if age < 90 days then try to load it. If loading fails retrieve data
        # from the NASA server .
        r = os.stat(cache_file)
        cache_file_date = date.fromtimestamp(r.st_mtime)
        age = (date.today() - cache_file_date).days
        if age < 90:
            status = self._load_cache_file()
            if status is not True:
                # Loading cache file failed!
                self._get_and_process_NASAPower(self.latitude, self.longitude)
        else:
            # Cache file is too old. Try loading new data from NASA
            try:
                self._get_and_process_NASAPower(self.latitude, self.longitude)
            except Exception, e:
                msg1 = "Failed updating weather data from NASA Power due to: %s" % e
                msg2 = "Falling back on outdated cache file '%s'" % cache_file
                self.logger.warning([msg1, msg2])


    def _get_and_process_NASAPower(self, latitude, longitude):
        """Handles the retrieval and processing of the NASA Power data
        """
        powerdata = self._query_NASAPower_server(latitude, longitude)

        # Store the informational header then parse variables
        self.description = self._parse_header(powerdata)
        self.elevation = self._parse_elevation(powerdata)
        recs = self._process_power_records(powerdata)

        # Start building the weather data containers
        self._make_WeatherDataContainers(recs)

        # dump contents to a cache file
        cache_filename = self._get_cache_filename(latitude, longitude)
        self._dump(cache_filename)

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
        req = urllib.urlopen(url)

        if req.getcode() != self.HTTP_OK:
            msg = "Failed retrieving POWER data from %s" % server
            self.logger.error(msg)
            raise PCSEError(msg)

        powerdata = req.readlines()
        req.close()

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
        """Constructs the filename used for cache files given latitude longitude
        """
        fname = "%s_LAT%04i_LON%04i.cache" % (self.__class__.__name__,
                                              int(latitude), int(longitude))
        cache_filename = os.path.join(global_settings.METEO_CACHE_DIR, fname)
        return cache_filename

    def _write_cache_file(self):
        """Writes the meteo data from NASA Power to a cache file.
        """
        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        try:
            self._dump(cache_filename)
        except (IOError), excp:
            msg = "Failed to write cache to file: %s" % cache_filename
            self.logger.warning(msg)

    def _load_cache_file(self):
        """Loads the data from the cache file. Return True if successful.
        """
        cache_filename = self._get_cache_filename(self.latitude, self.longitude)
        try:
            self._load(cache_filename)
            return True
        except (EnvironmentError, EOFError), e:
            msg = "Failed to load cache from file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
            return False


    def _make_WeatherDataContainers(self, recs):
        """Create a WeatherDataContainers from recs, compute ET and store the WDC's.
        """

        for rec in recs:
            wdc = WeatherDataContainer(LAT=self.latitude, LON=self.longitude,
                                       ELEV=self.elevation)
            wdc.DAY = rec["day"]
            for pcse_name, power_name, conv, unit in self.pcse_variables:
                value = conv(rec[power_name])
                wdc.add_variable(pcse_name, value, unit)

            # Reference evapotranspiration in mm/day
            try:
                (E0, ES0, ET0) = penman(wdc.DAY, wdc.LAT, wdc.ELEV, self.angstA,
                                        self.angstB, wdc.TMIN, wdc.TMAX, wdc.IRRAD,
                                        wdc.VAP, wdc.WIND)
            except ValueError, exc:
                msg = (("Failed to calculate reference ET values on %s" % wdc.DAY) +
                       "With input values:\n %s" % str(wdc))
                raise PCSEError(msg)

            # Add to wdc and convert to cm/day
            wdc.add_variable("E0", E0 / 10., "cm/day")
            wdc.add_variable("ES0", ES0 / 10., "cm/day")
            wdc.add_variable("ET0", ET0 / 10., "cm/day")

            # add wdc to dictionary for thisdate
            self._store_WeatherDataContainer(wdc, wdc.DAY)

    def _parse_elevation(self, powerdata):
        """Parse elevation out of the powerdata header"""
        for line in powerdata:
            if line.startswith("Elevation"):
                elev = int(line.split("=")[-1])
                return elev

    def _parse_header(self, powerdata):
        header = powerdata[1:7]
        header.append(powerdata[20])
        return [h.strip() for h in header]

    def _process_power_records(self, powerdata):
        """Process the meteorological records returned by NASA POWER
        """

        # First strip off the header by searching for '-END HEADER'
        is_header = True
        while is_header == True:
            line = powerdata.pop(0)
            if line.startswith("-END HEADER"):
                is_header = False

        # Now start parsing meteo records
        recs = []
        for line in powerdata:
            rec = self._parse_raw_power_record(line)
            if rec is not None:
                recs.append(rec)
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
            return None


