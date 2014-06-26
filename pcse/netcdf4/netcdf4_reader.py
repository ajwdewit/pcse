# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), May 2014
import os
import glob;
from netCDF4 import Dataset;
from datetime import date, datetime, time, timedelta;
import math;

from ..util import penman
from ..base_classes import WeatherDataProvider, WeatherDataContainer
from pcse.pcse.geo import AsciiGrid, GridEnvelope2D;
from ..exceptions import PCSEError
from ..settings import settings
import warnings

def ea_from_tdew(tdew):
    """ Arguments:
    tdew - dewpoint temperature [deg C]
    """
    if (tdew < -95.0 or tdew > 65.0):
        # Are these reasonable bounds?
        raise ValueError('tdew=%g is not in range -95 to +60 deg C' % tdew)
    
    # � la Priestley-Taylor and David Bolton
    tmp = (17.32491 * tdew) / (tdew + 238.102); 
    ea = 0.610588 * math.exp(tmp);
    # slope = 238.102*17.32491*VPS/(TMA+238.102)**2
    return ea; # in kPa

def vap_from_sh(sh, alt):
    """ Vapour pressure in hPa, from specific humidity 
    Arguments:
    sh  - specific humidity 
    alt - altitude above sea level
    """
    # Correct pressure for altitude
    p = 101.325 * math.pow(1 - 2.25577e-5 * alt, 5.2588)
    
    # Specific humidity is the mass of water vapour per mass of air mixture; therefore:
    mr = sh / (1 - sh); # mixing ratio, mass of water vapour per mass of dry air
    vap = p * mr / (0.622 + mr); # constant is based on ratio of mole weights for dry air and water
    return vap;

class NetcdfEnvelope2D(GridEnvelope2D):
    # Constants
    LON = 'lon';
    LAT = 'lat';
    TIME = 'time'
    
    # File names
    elevation_grid = r"../geodata/glob_elevation_resampled.asc";
    number_grid = r"../geodata/grid_50deg_gld.asc";
    
    def __init__(self, ds):
        pass;
    
    @staticmethod
    def getEnvelope(ds):
        # Retrieve the x and y ranges stored in the netCDF dataset
        LON = NetcdfEnvelope2D.LON;
        LAT = NetcdfEnvelope2D.LAT;
        _dims = ds.dimensions;
        _vars = ds.variables;
        x_range= NetcdfEnvelope2D._readRange(_vars[LON]);
        y_range = NetcdfEnvelope2D._readRange(_vars[LAT]);
        
        # Now create an envelope object
        dx = GridEnvelope2D._getStep(x_range[0], x_range[1], len(_dims[LON]));
        dy = GridEnvelope2D._getStep(y_range[0], y_range[1], len(_dims[LAT]));
        xll = min(_vars[LON]) - 0.5*dx;
        yll = min(_vars[LAT]) - 0.5*dy;
        nvlp = GridEnvelope2D(len(_dims[LON]), len(_dims[LAT]), xll, yll, dx, dy);
        
        # In some netCDF files the longitudes and latitudes are stored differently
        if _vars[LON][0] > _vars[LON][-1]: nvlp.xcoords_sort = 'DESC';
        if _vars[LAT][0] < _vars[LAT][-1]: nvlp.ycoords_sort = 'ASC';
        return nvlp;
    
    @staticmethod
    def _readRange(varxy):
        # Argument is either x or y
        minxy = min(varxy);
        maxxy = max(varxy);
        return [minxy, maxxy];
    
    @staticmethod
    def _readDimensions(ncdimensions):
        LON = NetcdfEnvelope2D.LON;
        LAT = NetcdfEnvelope2D.LAT;
        nrows = len(ncdimensions[LAT]);
        ncols = len(ncdimensions[LON]);
        return nrows, ncols;
    

class NetcdfWeatherDataProvider(WeatherDataProvider, NetcdfEnvelope2D):
    """WeatherDataProvider for using netcdf4 files with PCSE
    
    :param longitude: longitude to request weather data for
    :param latitude: latitude to request weather data for
    
    Weather data can efficiently be delivered in the form of files in the NETCDF4
    format. In general such data pertain to gridded weather, meaning that they are
    interpolated for a grid with regular intervals in space and time.
    
    We assume that the grid conforms to one of the WGS standards with longitude 
    representing the west to east direction and latitudes the north to south
    direction.
    
    This class can be developed further into a base class for use with different
    datasets. Such datasets have different names, start year, end year as well as
    different variables for representing the weather.
    
    """

    # Define some lambda functions to take care of unit conversions.
    W_to_J_day = lambda x: x * 86400;
    KtoC = lambda x : x - 273.15
    no_conv = lambda x: x
    Kg_M2_Sec_to_cm_day = lambda x: 86400 * x/10.; 
    #rh_to_hpa = lambda rh, t: rh * ea_from_tdew(t) / 10.;
    #sh_to_hpa = lambda sh, alt: vap_from_sh(sh, alt);
    
    # Variable names in dataset - append "hur" or "hus" in subclass
    netcdf_variables = ["pr", "rsds", "tas", "tasmin","tasmax", "wind"];
    
    # Mapping PCSE name to power name, conversion factor and unit of weather variables
    # Add optional tuples in the subclass - e.g. ("VAP", "hus", sh_to_hpa, "hPa")
    pcse_variables = [("IRRAD", "rsds", W_to_J_day, "J/m2/day"),
                      ("TMIN", "tasmin", KtoC, "Celsius"), 
                      ("TMAX", "tasmax", KtoC, "Celsius"),
                      ("TEMP", "tas", KtoC, "Celsius"),
                      ("WIND", "wind", no_conv, "m/sec"),
                      ("RAIN", "pr", Kg_M2_Sec_to_cm_day, "cm/day")];
 
    # Dictionary to store the relevant datasets:
    Netcdf4_files = {};
    available_years = [];

    # other constants
    angstA = 0.25
    angstB = 0.5

    # Pls note that we use lon and then lat, just because x is usu. mentioned before y
    def __init__(self, fname, longitude, latitude, fpath=None, force_update=False):
        WeatherDataProvider.__init__(self);
        
        # Construct search path
        search_path = self._construct_search_path(fname, fpath);
        
        # Of the required files, find out which are available 
        self.Netcdf4_files, self.available_years = self._get_Netcdf4_files(fname, search_path);

        # Access the first file and calculate the dimensions etc.
        key1 = self.Netcdf4_files.keys()[0]; # first key
        ds = self.Netcdf4_files[key1];
        nvlp = self.getEnvelope(ds);
        self.longitude, self.latitude = nvlp.getNearestCenterPoint(longitude, latitude);
        
        # Check for existence of a cache file
        cache_file = self._find_cache_file(self.longitude, self.latitude);
        if cache_file is None or force_update is True:
            msg = "No cache file or forced update, retrieving data from disk."
            self.logger.debug(msg)
            # No cache file, we really have to get the data from disk
            self._get_and_process_Netcdf4(fname, self.longitude, self.latitude, nvlp);
            return;       

        # Get age of cache file, if any of the Netcdf4 files is younger then try 
        # to load it. If loading fails, still try to retrieve data from the Netcdf4 files
        rc = os.stat(cache_file);
        cache_file_date = datetime.fromtimestamp(rc.st_mtime);
        agediff = -36525; # about 100 years ago
        
        # Do we have to re-establish connection to the datasets?
        # self.Netcdf4_files, self.available_years = self._get_Netcdf4_files(fname, search_path);
        for key in self.Netcdf4_files:
            # Check the file dates against the date of the cache file
            ds = self.Netcdf4_files[key];
            rnc4 = os.stat(ds.filepath());
            nc4_file_date = datetime.fromtimestamp(rnc4.st_mtime);
            agediff = max(agediff, (nc4_file_date - cache_file_date).days);
            
        if agediff < 0:
            # Cache file
            msg = "Start loading weather data from cache file: %s" % cache_file
            self.logger.debug(msg);
            status = self._load_cache_file();
            if status is not True:
                msg = "Loading cache file failed, reloading data from Netcdf4 files."
                self.logger.debug(msg)
                # Loading cache file failed!
                self._get_and_process_Netcdf4(self.longitude, self.latitude);
        else:
            # Cache file is too old. Try loading new data from file
            try:
                msg = "Cache file older then the Netcdf4 files, reloading data."
                self.logger.debug(msg)
                self._get_and_process_Netcdf4(self.longitude, self.latitude);
            except:
                msg = ("Reloading data from Netcdf4 files failed, reverting to (outdated) " +
                       "cache file")
                self.logger.debug(msg)
                status = self._load_cache_file()
                if status is not True:
                    msg = "Outdated cache file failed loading."
                    raise PCSEError(msg);

    
    def _get_Netcdf4_files(self, fname, search_path):
        """Find Netcdf4 files on given path with given name
        Also sorts the list, checks for missing years and sets self.lastyear
        Assume this pattern for the ncfile names:
        [variable_name]_[dataset_name.lower()]_[19??]-20??.[nc4]
        """
        if not search_path.endswith(os.path.sep): search_path += os.path.sep;
        tmpnc4_files = sorted(glob.glob(search_path + "*_" + fname.lower() + "_19??-20??.nc4"));
        if len(tmpnc4_files) == 0:
            msg = "No Netcdf4 files found when searching at %s"
            raise PCSEError(msg % search_path);
        
        mainresult = {};
        for ncfile in tmpnc4_files:
            # Check that we need this ncfile in the first place
            fn = os.path.basename(ncfile);
            pos = fn.find('_'); # it is assumed that the first part indicates the variable
            varname = fn[0:pos];
            if (varname in self.netcdf_variables):
                # Ok, the ncfile contains data wrt. one of the relevant variables
                ds = Dataset(ncfile, 'r');
                mainresult[varname] = ds; 
            
        # Retrieve the years from the ncfile names, taking into account that there
        # may be an alternative end_year
        start_year = 1900;
        for varname in mainresult:
            fn = os.path.basename(mainresult[varname].filepath());
            pos1 = str(fn).rfind('_');
            pos2 = str(fn).rfind('-');
            tmpyr = int(fn[pos1+1:pos2]);
            start_year = max(start_year, tmpyr);
        
        end_year = 2100;
        for varname in mainresult:
            fn = os.path.basename(mainresult[varname].filepath());
            pos1 = str(fn).rfind('-');
            pos2 = str(fn).rfind('.');
            tmpyr = int(fn[pos1+1:pos2]);
            end_year = min(end_year, tmpyr);
        
        available_years = range(start_year, end_year + 1);
        self.firstyear = start_year;
        self.lastyear = end_year;
        
        return mainresult, available_years;
    
    def _get_and_process_Netcdf4(self, fname, longitude, latitude, nvlp1):      
        # First check that the files do cover the same extent
        dateref = date(1860, 1, 1);
        timeref = datetime.combine(dateref, time(0,0,0));
        
        # Establish the extent for the first file
        key1 = self.Netcdf4_files.keys()[0]; # first key
        ds = self.Netcdf4_files[key1];
                
        # Compare the extent of the given envelope with that of the other one
        for key in self.Netcdf4_files:
            # It was already checked that the file contains data wrt. a relevant variable
            if (key == key1): continue; # given envelope based on content of first file
            ds = self.Netcdf4_files[key];
            nvlp = self.getEnvelope(ds);

            # Now check the extent and sorting of coordinates; if it's ok, then continue
            if not nvlp.hasSameExtent(nvlp1):
                raise PCSEError("Netcdf4 files do not cover the same extent");
            if not nvlp.compareSorting(nvlp1):
                raise PCSEError("Netcdf4 files do not have their coordinates sorted in the same way");

        # If we reach here, we can assume that all file shave the same extent
        self.elevation = self.get_elevation(longitude, latitude);
        self.description = "Meteo data from Netcdf4 files with label " + fname;
        
        # Prepare to check the dates for which data are available in the various files
        start_day = date(self.available_years[0], 1, 1);
        end_day = date(self.available_years[-1], 12, 31);
        
        # Get hold of all the data relevant for the given location, for the available years
        # Lookup lat-lon in the netCDF and check that latitudes and longitudes are the same! 
        ds.variables[self.LAT]
        k, i = nvlp1.getColAndRowIndex(self.longitude, self.latitude);
        latitudes = ds.variables[self.LAT];
        longitudes = ds.variables[self.LON];
        abs_diff = abs(latitudes[i] - self.latitude);
        assert abs_diff < 0.01, "Latitudes not equal: " + str(latitudes[i]) + " " + str(self.latitude);
        abs_diff = abs(longitudes[k] - self.longitude)
        assert abs_diff < 0.01, "Longitudes not equal!: " + str(longitudes[i]) + " " + str(self.longitude);

        dataslices = {}
        #for key in self.Netcdf4_files:
        #    ds = None;
        while self.Netcdf4_files:
            try:
                # Get hold of file and file name
                #ds = self.Netcdf4_files[key];
                key, ds = self.Netcdf4_files.popitem();
                
                # Check the time reference
                times = ds.variables[self.TIME];
                
                if (times.units != u'days since ' + str(timeref)):
                    fn = os.path.basename(ds.filepath());
                    raise PCSEError("Time in file " + fn + " not expressed as days since " + str(timeref));
                
                # Get the first and last day for which there are data in the file and check
                if (timeref + timedelta(times[0]) - datetime.combine(start_day, time(0,0,0))).days > 0: 
                    msg = "File " + fn + " does not have data for days as early as " + str(start_day.date());
                    raise PCSEError(msg);
                if (timeref + timedelta(times[-1]) - datetime.combine(end_day, time(0,0,0))).days < 0:
                    msg = "File " + fn + " does not have data for all days until " + str(end_day.date());
                    raise PCSEError(msg);
                
                # If we reach here, then get slices of the data - assume data are available for each day!
                # We'll copy the data year by year to avoid the ValueError "array is too big"
                self.logger.info("Shape of the data for variable " + key + ": " + str(ds.variables[key].shape));
                numdays = ((timeref + timedelta(times[-1])) - (timeref + timedelta(times[0]))).days;
                # t1 = (timeref + timedelta(times[0]) - start_day).days; 
                # tn = numdays - (timeref + timedelta(times[-1]) - end_day).days - 1;
                dataslice = ds.variables[key][:, i, k];
                dataslices[key] = dataslice;  
            except Exception as e:
                fn = os.path.basename(ds.filepath());
                raise PCSEError("An error occurred while reading file " + fn + " (" + str(e) + ")");
            finally:
                ds.close();

        # Estimate Angstrom coefficients
        self.AngstA, self.AngstB = self.getAngstromCoeff(self.latitude);

        # Prepare to loop over all the days in the dataset
        for day in range(0, numdays):
            t = {"LAT": self.latitude, "LON": self.longitude, "ELEV": self.elevation}
            t["DAY"] = start_day + timedelta(day);
            for pcse_name, netcdf_name, conv, units in self.pcse_variables:
                dataslice = dataslices[netcdf_name];
                if pcse_name != "VAP":
                    value = conv(dataslice[day]);  
                else:
                    if (netcdf_name =="hur"):
                        value = dataslice[day] * ea_from_tdew(t["TMIN"]) / 10.; # rh_to_hpa
                    else:
                        # convert specific humidity "hus" to vapour pressure
                        value = vap_from_sh(dataslice[day], t["ELEV"]); # sh_to_hpa
                t[pcse_name] = value;

            # Reference evapotranspiration in mm/day
            try:
                (E0,ES0,ET0) = penman(t["DAY"], t["LAT"], t["ELEV"], self.angstA,
                                      self.angstB, t["TMIN"], t["TMAX"], t["IRRAD"],
                                      t["VAP"], t["WIND"])
            except ValueError as e:
                msg = (("Failed to calculate reference ET values on %s. " % t["DAY"]) +
                       ("With input values:\n %s.\n" % str(t)) +
                       ("Due to error: %s" % e))
                raise PCSEError(msg)

            # update record with ET values value convert to cm/day
            t.update({"E0": E0/10., "ES0": ES0/10., "ET0": ET0/10.})

            # Build weather data container from dict 't'
            wdc = WeatherDataContainer(**t)

            # add wdc to dictionary for this date
            self._store_WeatherDataContainer(wdc, wdc.DAY);
        self._write_cache_file();
    
    def _find_cache_file(self, longitude, latitude):
        """Try to find a cache file for given latitude/longitude.
        Returns None if the cache file does not exist, else it returns the full path
        to the cache file.
        """
        cache_filename = self._get_cache_filename(longitude, latitude)
        if os.path.exists(cache_filename):
            return cache_filename
        else:
            return None;
    
    def _get_cache_filename(self, longitude, latitude):
        """Constructs the filename used for cache files given latitude and longitude
        The latitude and longitude is coded into the filename - no truncating.So the
        cache filename for a point with lat/lon 52.75/-124.75 will be:
        NetcdfWeatherDataProvider_LAT05275_LON-12475.cache
        """
        clsName = self.__class__.__name__;
        fname = "%s_LAT%005i_LON%005i.cache" % (clsName, int(latitude*100), int(longitude*100))
        cache_filename = os.path.join(settings.METEO_CACHE_DIR, fname);
        return cache_filename;
    
    def _write_cache_file(self):
        """Write the data loaded from the Netcdf files to a binary file using cPickle
        """
        cache_filename = self._get_cache_filename(self.longitude, self.latitude)
        try:
            self._dump(cache_filename)
        except (IOError, EnvironmentError), e:
            msg = "Failed to write cache to file '%s' due to: %s" % (cache_filename, e)
            self.logger.warning(msg)
    
    def _construct_search_path(self, fname, fpath):
        """Construct the path where to look for files"""
        if fpath is None:
            # assume NC4 files in current folder
            p = os.path.join(os.getcwd(), fname);
        elif os.path.isabs(fpath):
            # absolute path specified
            p = os.path.join(fpath, fname);
        else:
            # assume path relative to current folder
            p = os.path.join(os.getcwd(), fpath, fname);

        return os.path.normpath(p);
    
    def _load_cache_file(self): 
        """Load the weather data from a binary file using cPickle.

        Also checks if any of the Netcdf4 files have modification/creation date more recent then the cache_file.
        In that case reload the weather data from the original Netcdf files.
        
        Returns True if loading succeeded, False otherwise
        """
        # If no cache_file defined return False directly
        cache_filename = self._get_cache_filename(self.longitude, self.latitude);
        if not os.path.exists(cache_filename):
            return False;
        
        # Retrieve file dates
        nc_dates = [];
        for key in self.Netcdf4_files:
            ds = self.Netcdf4_files[key];
            nc_file = ds.filepath();
            r = os.stat(nc_file);
            nc_dates.append(r.st_mtime);

        # retrieve modification dates of cache file
        cache_date = os.stat(cache_filename).st_mtime
        if any([ncd > cache_date for ncd in nc_dates]):
            try:
                os.remove(cache_filename)
            except OSError, exc:
                msg = "Failed to remove cache file '%s' due to: %s" % (cache_filename, exc)
                warnings.warn(msg)
            return False
        else:
            # Else load data from cache file and store internally
            try:
                self._load(cache_filename)
                return True
            except:
                msg = "Cache file failed loading! Try to delete cache file: %s"
                self.logger.warn(msg, cache_filename)
                return False;
            
    def getAngstromCoeff(self, deglat):
        # Constants from pcse.util
        MIN_A = 0.1
        MAX_A = 0.4
        MIN_B = 0.3
        MAX_B = 0.7;
        
        # See van der Drift en van Diepen, 1992
        result = [0.25, 0.5]; # default values
        try:
            # Initialise A and B
            A = result[0];
            B = result[1];
            
            # THe folowing is based on method check_angstromAB found in pcse.util
            if abs(deglat) < 45.0:
                # For tropical regions, the A becomes too high and the B too low
                A = min(0.4885 - 0.0052 * abs(float(deglat)), MAX_A);
                B = max(0.1563 + 0.0074 * abs(float(deglat)), MIN_B);
            else:
                # For temperate regions, the A becomes too low and the B too high
                A = max(0.4885 - 0.0052 * abs(float(deglat)), MIN_A);
                B = min(0.1563 + 0.0074 * abs(float(deglat)), MAX_B);

            # Assign the result - for no latitude the sum of A and B will be too large
            result = [A, B];
        finally:
            return result;
    
    def _get_value_from_grid(self, longitude, latitude, fpath):
        # Open the file. Get right row and column. Elevations are linked to the cell centres
        r = AsciiGrid(fpath, "i");
        if not r.open('r'): raise Exception("Unable to open input file " + r.name);
        k, i = r.getColAndRowIndex(longitude, latitude);
        if (i == r.nrows): i = i - 1;
        
        # Now get hold of the right row, read the wanted value and close
        for _ in range(0, i): r.next(False);
        line = r.next(); # this line is split, unlike the previous ones
        r.close();
        return line[int(k)];    
    
    def get_elevation(self, longitude, latitude):
        # Find out where the elevation grid might be located
        fname = self.elevation_grid;
        fullpath = self._getFullPath(fname);
        return self._get_value_from_grid(longitude, latitude, fullpath);
    
    def get_grid_no(self, longitude, latitude):
        # Find out where the number grid might be located
        fname = self.number_grid;
        fullpath = self._getFullPath(fname);
        return self._get_value_from_grid(longitude, latitude, fullpath);
    
    def _getFullPath(self, fname):
        key1 = self.Netcdf4_files.keys()[0]; # first key
        ds = self.Netcdf4_files[key1];
        path = os.path.dirname(ds.filepath());
        return os.path.join(path, fname);
    
    def close(self):
        # finally close the files
        for key in self.Netcdf4_files:
            try:
                # Get hold of file and file name
                ds = self.Netcdf4_files[key];
                ds.close();
            finally:
                del(ds);


    