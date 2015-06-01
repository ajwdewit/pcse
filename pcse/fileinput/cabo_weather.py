# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import os, sys
import glob
import calendar
import numpy as np
import datetime as dt
import warnings

from ..base_classes import WeatherDataContainer, WeatherDataProvider
from ..util import reference_ET, angstrom, check_angstromAB
from ..exceptions import PCSEError

class CABOWeatherDataProvider(WeatherDataProvider):
    """Reader for CABO weather files.
    
    :param fname: root name of CABO weather files to read
    :param fpath: path where to find files, can be absolute or relative.
    :keyword ETmodel: "PM"|"P" for selecting penman-monteith or Penman
        method for reference evapotranspiration. Defaults to "PM".
    :keyword distance: maximum interpolation distance for meteorological
        variables, defaults to 1 day.
    :returns: callable like object with meteo records keyed on date.
    
    The Wageningen crop models that are written in FORTRAN or FST often use
    the CABO weather system (http://edepot.wur.nl/43010) for storing and
    reading weather data. This class implements a reader for the CABO weather
    files and also implements additional features like interpolation of
    weather data in case of missing values, conversion of sunshine duration
    to global radiation estimates and calculation the reference
    evapotranspiration values for water, soil and plants (E0, ES0, ET0)
    using the Penman approach.
    
    A difference with the old CABOWE system is that the python implementation
    will read and store all files (e.g. years) available for a certain
    station instead of loading a new file when crossing a year boundary.
    
    .. note::
        some conversions are done by the CABOWeaterDataProvider from
        the units in the CABO weather file for compatibility with WOFOST:
        
        - vapour pressure from kPa to hPa
        - radiation from kJ/m2/day to J/m2/day
        - rain from mm/day to cm/day.
        - all evaporation/transpiration rates are also returned in cm/day.
    
    *Example*
    
    The file 'nl1.003' provides weather data for the year 2003 for the
    station in Wageningen and can be found in the cabowe/ folder of the
    WOFOST model distribution. This file can be read using::
    
        >>> weather_data = CABOWeatherDataProvider('nl1', fpath="./meteo/cabowe")
        >>> print weather_data(datetime.date(2003,7,26))
        Weather data for 2003-07-26 (DAY)
        IRRAD:  12701000.00  J/m2/day
         TMIN:        15.90   Celsius
         TMAX:        23.00   Celsius
          VAP:        16.50       hPa
         WIND:         3.00     m/sec
         RAIN:         0.12    cm/day
           E0:         0.36    cm/day
          ES0:         0.32    cm/day
          ET0:         0.31    cm/day
        Latitude  (LAT):    51.97 degr.
        Longitude (LON):     5.67 degr.
        Elevation (ELEV):    7.0 m.

    Alternatively the date in the print command above can be specified as
    string with format YYYYMMDD or YYYYDDD.
    """
    # Order, name and conversion factor of weather variables
    variables = [("IRRAD",1e3, "J/m2/day"),("TMIN",1,"Celsius"),
                 ("TMAX",1,"Celsius"),("VAP",10,"hPa"),
                 ("WIND",1,"m/sec"), ("RAIN",0.1,"cm/day")]
    # Radiation in Kj/m2/day or Sunshine duration
    has_sunshine = False
    # Status line and observation no data value
    status_no_data = -999.
    weather_no_data = -99.
    #  start and end year
    firstyear = None
    lastyear = None
    # First date when CABO files start
    first_date = None
    # temporary array for storing data
    potential_records = None
    tmp_data = None
    
    def __init__(self, fname, fpath=None, ETmodel="PM", distance=1):
        WeatherDataProvider.__init__(self)

        self.ETmodel = ETmodel

        # Construct search path
        search_path = self._construct_search_path(fname, fpath)
        # find available files
        CABOWE_files, available_years, cache_file = self._find_CABOWEfiles(search_path)

        # If no cache file can be found, then start loading the CABOWE files
        if not self._load_cache_file(cache_file, CABOWE_files):
            self.tmp_data = self._calc_arraysize()

            # Run through files, read header and location parameters.
            # Then read meteo data into tmp_data array
            prev_cb_file = None
            for yr, cb_file in zip(available_years, CABOWE_files):
                header, loc_par, records = self._read_file(cb_file)
                # header info is taken from the first CABOWE file
                if self.description is None:
                    self.description = header
                self._set_location_parameters(loc_par, cb_file, prev_cb_file)

                for rec in records:
                    if len(rec.strip()) == 0:
                        continue
                    if rec.startswith("-999"):
                        # Status line without values
                        continue
                    self._proc_weather_record(rec, yr)
                prev_cb_file = cb_file

            # convert sunshine duration to global radiation
            self._check_angstrom()
            # Run interpolation, given interpolation distance (default is 1 day)
            self._interpolate_timeseries(distance)
            self._make_WeatherDataContainers()

            # Write data to binary cache file
            self._write_cache_file(search_path)

            # Delete array for tmp storage
            delattr(self, "tmp_data")

    def _load_cache_file(self, cache_file, CABOWE_files):
        """Load the weather data from a binary file using cPickle.

        Also checks if any of the CABOWE files have modification/creation date more recent then the cache_file.
        In that case reload the weather data from the original CABOWE files.
        
        Returns True if loading succeeded, False otherwise
        """
        # If no cache_file defined return False directly
        if cache_file is None:
            return False

        # if date of any CABOWE files > cache file: discard cache file and return False to reload from original files.
        # retrieved last modification dates of CABOWE files
        cb_dates = []
        for cb_file in CABOWE_files:
            r = os.stat(cb_file)
            cb_dates.append(r.st_mtime)
        # retrieve modification dates of cache file
        cache_date = os.stat(cache_file).st_mtime
        if any([cbd > cache_date for cbd in cb_dates]):
            try:
                os.remove(cache_file)
            except OSError as exc:
                msg = "Failed to remove cache file '%s' due to: %s" % (cache_file, exc)
                warnings.warn(msg)
            return False
        else:
            # Else load data from cache file and store internally
            try:
                self._load(cache_file)
                return True
            except Exception as e:
                msg = "Cache file failed loading! Try to delete cache file: %s"
                self.logger.warn(msg, cache_file)
                return False

    def _write_cache_file(self, search_path):
        """Write the data loaded from the CABOWE files to a binary file using cPickle
        """
        cache_fname = search_path + ".cache"
        self._dump(cache_fname)

    def _construct_search_path(self, fname, fpath):
        """Construct the path where to look for files"""
        if fpath is None:
            # assume CABOWE files in current folder
            p = os.path.join(os.getcwd(), fname)
        elif os.path.isabs(fpath):
            # absolute path specified
            p = os.path.join(fpath, fname)
        else:
            # assume path relative to current folder
            p = os.path.join(os.getcwd(), fpath, fname)

        return os.path.normpath(p)

    def _proc_weather_record(self, rec, fileyr):
        """Processes record and inserts values into correct place in array"""
        values = rec.split()
        try:
            year = int(values[1])
            doy  = int(values[2])
            weather_obs = np.array(values[3:], dtype=np.float64)
        except (ValueError,IndexError) as exc:
            msg = ("Failed to parse line: %s" % rec)
            raise RuntimeError(msg)
        
        # Check if file contents are consistent with file year
        if year != fileyr:
            msg = "File with year %s contains record for year %s"
            raise PCSEError(msg % (fileyr, year))

        # Calculate position in tmp_array base on date since first date
        rec_date = dt.date(year, 1, 1) + dt.timedelta(days=(doy-1))
        arraypos = (rec_date - self.first_date).days
        
        # insert data at right position in array
        self.tmp_data[:, arraypos] = weather_obs
            
    def _interpolate_timeseries(self, distance=1):
        """Interpolates gaps using linear interpolation, except for rainfall.
        
        distance specifies the interpolation distance: distance=1 will only
        allow interpolation when the previous and the following day are
        available, distance=2 allows also when 2 days are missing, etc.
        Defaults to 1.
        """
        # Kernel for interpolation distance
        kernel = np.ones((1 + distance*2))

        # array for tracing missing values in the self.tmp_data
        has_data = np.ones_like(self.tmp_data)
        # find indivual missing observations which are equal to self.weather_no_data (e..g -99.),
        # set them to np.NaN
        index = np.where(self.tmp_data == self.weather_no_data)
        has_data[index] = 0
        self.tmp_data[index] = np.NaN
        # Find missing lines in the CABOWE files, these have not been inserted in tmp_data and
        # are therefore np.NaN (tmp_data was initialized with np.NaN)
        index = np.where(np.isnan(self.tmp_data) == True )
        has_data[index] = 0

        for i, (var, cf, unit) in enumerate(self.variables):
            if var == "RAIN":
                # No interpolation on rainfall data
                continue
            timeseries_hasdata = has_data[i,:].flatten()
            if timeseries_hasdata.sum() == timeseries_hasdata.size:
                # No missing values
                continue
            timeseries = self.tmp_data[i,:].flatten()
            r = np.convolve(timeseries_hasdata, kernel, mode='same')
            
            # Find positions for interpolation: hasdata==0 and >=2 neighbours
            # except the first and last record
            index = (timeseries_hasdata==0)*(r>=2)
            index[0]  = False
            index[-1] = False
            if True not in index:
                # No positions that can be interpolated
                continue
                
            # Determine positions and y values for interpolation (x, xp, yp)
            xrange = np.arange(self.potential_records, dtype=np.float)
            x  = xrange[index]
            xp = xrange[(timeseries_hasdata == 1)]
            yp = timeseries[(timeseries_hasdata == 1)]
            y_int = np.interp(x,xp,yp)

            # put interpolated values back into tmp_data
            self.tmp_data[i, x.astype(np.int)] = y_int
    

    def _make_WeatherDataContainers(self):
        """Converts the data in self.tmp_data into WeatherDataContainers which are stored in
        the class dictionary keyed on the date.
        
        Records that are incomplete (contain np.NaN values) are skipped.
        Moreover, if the radiation measurements are in sunshine duration, than
        the Angstrom equation is used to estimate global radiation. Finally,
        the evapotranspiration value are calculated for each complete record.
        """
        
        # Generate prototype weather data container
        #wdc_proto = self._build_WeatherDataContainer()
        
        for i in range(self.potential_records):
            rec = self.tmp_data[:, i]
            if True in np.isnan(rec):
                # Incomplete record: skip
                continue

            # Derive date from position in array
            thisdate = self.first_date + dt.timedelta(days=i)
            t = {"DAY": thisdate, "LAT": self.latitude, 
                 "LON": self.longitude, "ELEV": self.elevation}
            
            for obs, (name, cf, unit) in zip(rec, self.variables):
                if name == "IRRAD" and self.has_sunshine is True:
                    obs = angstrom(thisdate, self.latitude, obs, self.angstA, self.angstB)
                    # angstrom routine returns in J/m2/day, no conversion factor needed
                    t[name] = float(obs)
                else:
                    t[name] = float(obs)*cf
            
            # Reference evapotranspiration in mm/day
            try:
                E0, ES0, ET0 = reference_ET(thisdate, t["LAT"], t["ELEV"], t["TMIN"], t["TMAX"], t["IRRAD"],
                                            t["VAP"], t["WIND"], self.angstA, self.angstB, self.ETmodel)
            except ValueError as e:
                msg = (("Failed to calculate reference ET values on %s. " % thisdate) +
                       ("With input values:\n %s.\n" % str(t)) +
                       ("Due to error: %s" % e))
                raise PCSEError(msg)

            # update record with ET values value convert to cm/day
            t.update({"E0": E0/10., "ES0": ES0/10., "ET0": ET0/10.})

            # Build weather data container from dict 't'
            wdc = WeatherDataContainer(**t)

            # add wdc to dictionary for thisdate
            self._store_WeatherDataContainer(wdc, thisdate)
        
    def _calc_arraysize(self):
        """Returns array of NaNs with size based on min/max year of data."""
        self.potential_records = 0
        for yr in range(self.firstyear, self.lastyear+1):
            if calendar.isleap(yr):
                self.potential_records += 366
            else:
                self.potential_records += 365
        t_ar = np.empty((6,self.potential_records), dtype=np.float64)
        t_ar[:] = np.NaN
        
        return t_ar

    def _set_location_parameters(self, line, cb_file, prev_cb_file):
        """Parse, check and assign location parameters.
        """
        strvalues = line.split()
        if len(strvalues) != 5:
            msg = "Did not find 5 values on location parameter line of file %s"
            raise PCSEError(msg % file)
        
        parnames = ["longitude", "latitude", "elevation", "angstA", "angstB"]
        for parname, strvalue in zip(parnames, strvalues):
            try:
                fvalue = float(strvalue)
                current_value = getattr(self, parname)
                if current_value is None:
                    setattr(self, parname, fvalue)
                else:
                    if abs(current_value - fvalue) > 0.001:
                        raise AttributeError
            except ValueError as e:
                msg = "Failed to parse location parameter %s on file %s, value: %s"
                raise PCSEError(msg % (parname, cb_file, strvalue))
            except AttributeError as e:
                msg = "Inconsistent '%s' location parameter in file %s compared to file %s."
                raise PCSEError(msg % (parname, cb_file, prev_cb_file))
    
    def _check_angstrom(self):
        """Checks the Angstrom parameters for consistency.
        
        Also sets self.has_sunshine=True when both A and B > 0.
        """
        if self.angstA > 0 and self.angstB > 0:
            self.has_sunshine=True
            
        self.angstA = abs(self.angstA)
        self.angstB = abs(self.angstB)
        check_angstromAB(self.angstA, self.angstB)

    def _read_file(self, fname):
        
        with open(fname) as fp:
            lines = fp.readlines()
        header = []
        location_par = None
        records = []
        for line in lines:
            l = line.strip()
            if l.startswith("*"):
                header.append(l)
            else:
                if location_par is None:
                    location_par = l
                else:
                    records.append(l)
        return (header, location_par, records)
    
    def _find_CABOWEfiles(self, search_path):
        """Find CABOWEfiles on given path with given name.
        
        Also sorts the list, checks for missing years and sets self.firstyear/lastyear
        and self.first_date.
        """
        
        cachefile = search_path + ".cache"
        if not os.path.exists(cachefile):
            cachefile = None

        CABOWEfiles = sorted(glob.glob(search_path+".[0-9][0-9][0-9]"))
        if len(CABOWEfiles) == 0:
            path, tail = os.path.split(search_path+".???")
            msg = "No CABO Weather files found when searching for '%s' at %s"
            raise PCSEError(msg % (tail, path))

        available_years = []
        for Cfile in CABOWEfiles:
            path, ext = os.path.splitext(Cfile)
            ext = ext[1:]
            if ext.startswith("9"):
                weather_year = 1000 + int(ext)
            else:
                weather_year = 2000 + int(ext)
            available_years.append(weather_year)

        self.firstyear = min(available_years)
        self.lastyear = max(available_years)
        self.first_date = dt.date(self.firstyear, 1, 1)

        # Check if years are missing and if so write a warning to the log file
        all_years = set(range(self.firstyear, self.lastyear+1))
        diff = all_years.difference(set(available_years))
        if len(diff) > 0:
            msg = "No CABOWE files found for year(s): %s" % list(diff)
            warnings.warn(msg)

        return CABOWEfiles, available_years, cachefile
