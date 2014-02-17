import sys, os
import urllib
from datetime import date, timedelta

from .base_classes import WeatherDataProvider, WeatherDataContainer
from .util import wind10to2, ea_from_tdew, penman
from .exceptions import PCSEError

# Define some lambdas to take care of unit conversions.
MJ_to_J = lambda x: x*1e6
no_conv = lambda x: x
mm_to_cm = lambda x: x/10.
tdew_to_hpa = lambda x: ea_from_tdew(x)*10

class NASAPowerWeatherDataProvider(WeatherDataProvider):
    """
    """

    # Location parameters
    angstA = 0.25
    angstB = 0.5
    longitude = None
    latitude  = None
    elevation = None

    # Variable names in POWER data
    power_variables = ["toa_dwn","swv_dwn","lwv_dwn","T2M","T2MN","T2MX","RH2M","DFP2M","RAIN","WS10M"]

    # Mapping PCSE name to power name, conversion factor and unit of weather variables
    pcse_variables = [("IRRAD","swv_dwn", MJ_to_J, "J/m2/day"),
                      ("TMIN","T2MN",no_conv,"Celsius"),
                      ("TMAX","T2MX",no_conv,"Celsius"),
                      ("TEMP","T2M",no_conv,"Celsius"),
                      ("VAP","DFP2M",tdew_to_hpa,"hPa"),
                      ("WIND","WS10M",wind10to2,"m/sec"),
                      ("RAIN","RAIN",mm_to_cm,"cm/day")]
    # other constants
    HTTP_OK = 200

    def __init__(self, latitude=None, longitude=None, cache_path=None,
                 force_update=None):
        WeatherDataProvider.__init__(self)

        if latitude < -90 or latitude > 90:
            msg = "Latitude should be between -90 and 90 degrees."
            raise ValueError(msg)
        if longitude < -180 or longitude > 180:
            msg = "Longitude should be between -180 and 180 degrees."
            raise ValueError(msg)

        self.latitude = float(latitude)
        self.longitude = float(longitude)

        # build URL for retrieving data
        server = "power.larc.nasa.gov"
        t_url = ("http://{server}/cgi-bin/cgiwrap/solar/agro.cgi?" + 
                 "email=agroclim%40larc.nasa.gov&step=1&lat={lat}&lon={lon}" + 
                 "&ms=1&ds=1&ys=1984&me={month}&de={day}&ye={year}&p=toa_dwn&" +
                 "p=swv_dwn&p=lwv_dwn&p=T2M&p=T2MN&p=T2MX&p=RH2M&" + 
                 "p=DFP2M&p=RAIN&p=WS10M&submit=Submit")
        d = date.today()
        url = t_url.format(server=server, lat=latitude, lon=longitude,year=d.year,
                           month=d.month, day=d.day)
        req = urllib.urlopen(url)

        if req.getcode() != self.HTTP_OK:
            msg = "Failed retrieving POWER data from %s" % server
            self.logger.error(msg)
            raise PCSEError(msg)

        powerdata = req.readlines()
        req.close()

        # Store the informational header then parse variables
        self._header = self._get_header(powerdata)
        self.elevation = self._get_elevation(self._header)
        recs = self._process_power_variables(powerdata)

        # Start building the weather data containers
        self._make_WeatherDataContainers(recs)


    def _make_WeatherDataContainers(self, recs):
        """Builds a WeatherDataContainer and sets the location parameters.

        All records will deepcopy this model and add the weather variables on
        it.
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
                (E0,ES0,ET0) = penman(wdc.DAY, wdc.LAT, wdc.ELEV, self.angstA,
                                      self.angstB, wdc.TMIN, wdc.TMAX, wdc.IRRAD,
                                      wdc.VAP, wdc.WIND)
            except ValueError, exc:
                msg = (("Failed to calculate reference ET values on %s" % thisdate) +
                        "With input values:\n %s" % str(wdc))
                raise PCSEError(msg)

            # Add to wdc and convert to cm/day
            wdc.add_variable("E0", E0/10.,"cm/day")
            wdc.add_variable("ES0",ES0/10.,"cm/day")
            wdc.add_variable("ET0",ET0/10.,"cm/day")

            # add wdc to dictionary for thisdate
            self._store_WeatherDataContainer(wdc, wdc.DAY)



    def _get_elevation(self, header):
        "Parse elevation out of the system"
        for line in header:
            if line.startswith("Elevation"):
                elev = int(line.split("=")[-1])
                return elev

    def _get_header(self, powerdata):
        header = powerdata[1:7]
        header += powerdata[20]
        return header

    def _process_power_variables(self, powerdata):
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
            rec = self._parse_meteorecord(line)
            if rec is not None:
                recs.append(rec)
        return recs

    def _parse_meteorecord(self, line):
        """Parse each record and return as a dict, return None if parsing fails due to a missing variable
        """
        r = line.split()
        rec = {}
        year = int(r.pop(0))
        doy = int(r.pop(0))
        rec["day"] = date(year,1,1) + timedelta(days=(doy-1))
        try:
            for i, name in enumerate(self.power_variables):
                rec[name] = float(r[i])
            return rec
        except ValueError:
            return None



                                  
    
    
    
