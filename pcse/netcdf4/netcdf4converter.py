# -*- coding: utf-8 -*-
from ..geo import Netcdf4Envelope2D, Netcdf4Raster, Hdf5Raster, FloatingPointRaster;
from ..util import reference_ET;
from ..exceptions import PCSEError
import os, glob, array;
import logging;
from ..traitlets import Instance;
from math import pow, exp;
from datetime import date, datetime, time, timedelta;
from tables import Time32Col, Float32Col, IsDescription;

def ea_from_tdew(tdew):
    """ Arguments:
    tdew - dewpoint temperature [deg C]
    """
    if (tdew < -95.0 or tdew > 65.0):
        # Are these reasonable bounds?
        raise ValueError('tdew=%g is not in range -95 to +60 deg C' % tdew)
    
    # à la Priestley-Taylor and David Bolton
    tmp = (17.32491 * tdew) / (tdew + 238.102); 
    ea = 0.610588 * exp(tmp);  # slope = 238.102*17.32491*VPS/(TMA+238.102)**2
    return ea; # in kPa

def vap_from_sh(sh, alt):
    """ Vapour pressure in hPa, from specific humidity 
    Arguments:
    sh  - specific humidity 
    alt - altitude above sea level
    """
    # Correct pressure for altitude
    p = 101.325 * pow(1 - 2.25577e-5 * alt, 5.2588)
    
    # Specific humidity is the mass of water vapour per mass of air mixture; therefore:
    mr = sh / (1 - sh); # mixing ratio, mass of water vapour per mass of dry air
    vap = p * mr / (0.622 + mr); # constant is based on ratio of mole weights for dry air and water
    return vap;

# Use empty object as container
class Container(object):
    pass

class NetcdfWeatherDataConverter():
    # File names
    landmask_grid = r"../geodata/glob_landmask_resampled.flt";
    elevation_grid = r"../geodata/glob_elevation_resampled.flt";
    
    # Define some lambda functions to take care of unit conversions.
    W_to_J_day = lambda x: x * 86400;
    KtoC = lambda x : x - 273.15
    no_conv = lambda x: x
    Kg_M2_Sec_to_cm_day = lambda x: 86400 * x/10.; 
    
    # Variable names in dataset - append "hur" or "hus" in subclass
    netcdf_variables = ["pr", "rsds", "tas", "tasmin","tasmax", "wind"];
    
    pcse_variables = [("IRRAD", "rsds", W_to_J_day, "J/m2/day"),
                  ("TMIN", "tasmin", KtoC, "Celsius"), 
                  ("TMAX", "tasmax", KtoC, "Celsius"),
                  ("TEMP", "tas", KtoC, "Celsius"),
                  ("WIND", "wind", no_conv, "m/sec"),
                  ("RAIN", "pr", Kg_M2_Sec_to_cm_day, "cm/day")];
    
    # Dictionary to store the relevant datasets:
    Netcdf4rasters = {};
    available_years = [];
    __dataset_name = "dummy";
    __dataslices = None;
    logger = Instance(logging.Logger);
    
    def __init__(self, fname, fpath=None):
        try:
            # Initialise
            loggername = "%s.%s" % (self.__class__.__module__, self.__class__.__name__);
            self.logger = logging.getLogger(loggername);
            
            # Construct search path
            search_path = self._construct_search_path(fname, fpath);
            
            # Of the required files, find out which are available 
            self.Netcdf4rasters, self.available_years = self.__get_hdf5_filefname, search_path);
            
            # Process them
            self.__dataslices = {};
            self.__dataset_name = fname;
            self._get_and_process_Netcdf4(fname);
            
            self.logger.info("Conversion successful!");
        except Exception as e:
            self.logger.error(str(e));
        
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
        
    def __get_hdf5_fileself, fname, search_path):
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
        
        key_fp_dict = {};
        for ncfile in tmpnc4_files:
            # Check that we need this ncfile in the first place
            fn = os.path.basename(ncfile);
            pos = fn.find('_'); # it is assumed that the first part indicates the variable
            varname = fn[0:pos];
            if (varname in self.netcdf_variables):
                # Ok, the ncfile contains data wrt. one of the relevant variables
                nc4r = Netcdf4Raster(ncfile);
                if not nc4r.open('r'):
                    raise PCSEError("An error occurred while opening file " + nc4r.getFilePath());
                key_fp_dict[varname] = nc4r;
            
        # Retrieve the years from the ncfile names, taking into account that there
        # may be an alternative end_year
        start_year = 1900;
        for key in key_fp_dict:
            fn = os.path.basename(key_fp_dict[key].getFilePath());
            pos1 = str(fn).rfind('_');
            pos2 = str(fn).rfind('-');
            start_year = max(start_year, int(fn[pos1+1:pos2]));
        
        end_year = 2100;
        for key in key_fp_dict:
            fn = os.path.basename(key_fp_dict[key].getFilePath());
            pos1 = str(fn).rfind('-');
            pos2 = str(fn).rfind('.');
            end_year = min(end_year, int(fn[pos1+1:pos2]));
        
        available_years = range(start_year, end_year + 1);
        self.firstyear = start_year;
        self.lastyear = end_year;
        
        return key_fp_dict, available_years;

    
    def __get_and_process_hdf5self, fname):      
        # First check that the files do cover the same period and extent
        dateref = date(1860, 1, 1);
        timeref = datetime.combine(dateref, time(0,0,0));
        
        # Establish the extent for the first file
        key1 = self.Netcdf4rasters.keys()[0]; # first key
        nc4r1 = self.Netcdf4rasters[key1];
        nvlp1 = nc4r1.getEnvelope();
                
        # Carry out some checks on the various netCDF4 rasters
        self._check_Netcdf4rasters(timeref, nvlp1)

        # Get a landmask and elevation data
        fullpath = self._getFullPath(self.landmask_grid);
        fpr_lm= FloatingPointRaster(fullpath, 'i');
        if not fpr_lm.open('r'): raise Exception("Unable to open FLT file with landmask!");
        fullpath = self._getFullPath(self.elevation_grid);
        fpr_elev = FloatingPointRaster(fullpath, 'i');
        if not fpr_elev.open('r'): raise Exception("Unable to open FLT file with elevation data!");
        
        h5r = None;
        try:
            # Open the output file
            n = nvlp1;
            datadir = os.path.dirname(nc4r1.getFilePath());
            fpath = os.path.join(datadir, fname + r'_data.h5');
            h5r = Hdf5Raster(fpath);
            
            # Typically the hdf5vars are: [day. tmax, tmin, temp, rain, irrad, wind, e0, es0, et0]
            hdf5vars = getOrderedKeys(ncdf_weather.columns); 
            _units = ['date', 'Celsius', 'Celsius', 'Celsius', 'cm/day', "J/m2/day" 'm/s', 'cm/day', 'cm/day', 'cm/day'];
            h5r.open('w', n.ncols, n.nrows, n.xll, n.yll, n.dx, nc4r1.nodatavalue, 
                "ncdf_weather", variables=hdf5vars, units=_units);
            del n;
    
            # If we reach here, we can assume that all files apply to the same extent and period
            # Loop over the rows (latitudes) as a preparation for writing the data to a HDF5 file,
            # i.e. row by row! i is the row index
            for i in range(0, nvlp1.nrows):
                # Retrieve a data slice for the current row; get input for the next row
                rawline = fpr_lm.next();
                landmask = array.array('f', rawline); 
                rawline = fpr_elev.next();
                elevations = array.array('f', rawline); 
                
                # Initialise a structure for the output data with for each pixel in this row a place
                latitude = nc4r1.getVariables(Netcdf4Envelope2D.LAT)[i];
                self.logger.info("About to convert data for the next row, pertaining to latitude %s" % latitude);   
                data = [None] * nvlp1.ncols;                

                # Check whether there is land in this row in the first place
                if not self._any_land(landmask):
                    # Do as little as possible as far as this row is concerned
                    for key in self.Netcdf4rasters:
                        nc4r = self.Netcdf4rasters[key];
                        nc4r.next(False); 
                        self.__dataslices[key] = None;
                    continue;

                # Leave the netCDF4 specific stuff hidden by using class netcdf4raster 
                for key in self.Netcdf4rasters:
                    try:
                        # Get hold of the variable name and the corresponding raster, then get a slice of the data
                        nc4r = self.Netcdf4rasters[key];
                        self.__dataslices[key] = nc4r.next(); 
                    except Exception as e:
                        fn = os.path.basename(nc4r.getFilePath());
                        raise PCSEError("An error occurred while reading file " + fn + " (" + str(e) + ")");
                        
                # If we reach here, then get slices of the data - assume data are available for each day!
                self.logger.info("Shape of the data: " + str(nc4r1.getVariables(key1).shape));
                times = nc4r1.getVariables(Netcdf4Envelope2D.TIME);
                numdays = ((timeref + timedelta(times[-1])) - (timeref + timedelta(times[0]))).days;
            
                # Within the loop over the rows, we need a loop over the columns (longitudes)
                # Use a land mask and leave places in the data None for those pixels; k is the column index
                for k in range(0, nvlp1.ncols):
                    # Check that this pixel represents a land surface
                    if landmask[k] != 1: continue;
                    self.check_location(nc4r1, k, i);
                    recs = self._process_pixel(k, i, numdays, elevations[k]);

                    # Report how many days are missing
                    if len(recs) < numdays: 
                        longitude = nc4r1.getVariables(Netcdf4Envelope2D.LON)[k];
                        self.logger.warn("Less than expected number of records obtained for longitude %s" % longitude);
            
                    # First try to write 1 pixel
                    h5r.write(k, recs, ncdf_weather);
                    h5r.flush();
                
                    # Add the records to the array data, but if all are missing then don't
                    if len(recs) > 0: data[k] = recs;
                
                # Output for current row is complete
                # h5r.writenext(data, ncdf_weather);
                # h5r.flush();
        except Exception, e:
            self.logger.error(str(e));
        finally:
            if h5r != None: h5r.close();
    
    def _any_land(self, landmask):
        for k in range(0, len(landmask)):
            if landmask[k] == 1: return True;
        return False;
    
    def _process_pixel(self, colIndex, rowIndex, numdays, elevation):
        # Initialise
        k = colIndex;
        i = rowIndex;
        ts = None;
        
        try:
            # Estimate Angstrom coefficients
            key1 = self.Netcdf4rasters.keys()[0]; # first key
            nvlp1 = self.Netcdf4rasters[key1].getEnvelope();
            longitude, latitude = nvlp1.getXandYfromIndices(k, i);
            AngstA, AngstB = self._getAngstromCoeff(latitude);
    
            # Typically the hdf5vars are: [day. tmax, tmin, temp, rain, irrad, wind, e0, es0, et0]
            hdf5vars = getOrderedKeys(ncdf_weather.columns); 
    
            # Prepare to loop over all the days in the dataset; use class ncdf_weather for storing the data 
            ts = Container()
            start_day = datetime(self.available_years[0], 1, 1).date();
            for day in range(0, numdays):
                # Initialise record for this day
                curdate = start_day + timedelta(days=day);
                rec = [curdate.toordinal(), None, None, None, None, None, None, None, None, None, None];                    
                
                # Retrieve relevant data from the various data slices
                empty_slice_found = False;
                for pcse_name, netcdf_name, conv, units in self.pcse_variables:
                    if self.__dataslices[netcdf_name] != None:
                        timeseries = self.__dataslices[netcdf_name][:, k]; 
                        if str(timeseries[day]) != '--':                    
                            if pcse_name != "VAP":
                                value = conv(timeseries[day]);  
                            else:
                                if (netcdf_name =="hur"):
                                    pos = hdf5vars.index("tmin");
                                    value = timeseries[day] * ea_from_tdew(rec[pos]) / 10.; # rh_to_hpa
                                else:
                                    # convert specific humidity "hus" to vapour pressure
                                    value = vap_from_sh(timeseries[day], elevation); # sh_to_hpa
                            pos = hdf5vars.index(pcse_name.lower());
                            rec[pos] = value;
                    else: 
                        empty_slice_found = True;
                
                if not empty_slice_found and self._record_is_ok(rec, 3):
                    try:
                        # Reference evapotranspiration in mm/day
                        (E0,ES0,ET0) = reference_ET(curdate, latitude, elevation,
                            rec[hdf5vars.index("tmin")], rec[hdf5vars.index("tmax")],
                            rec[hdf5vars.index("irrad")], rec[hdf5vars.index("vap")],
                            rec[hdf5vars.index("wind")], AngstA, AngstB, "PM");
                    except ValueError as e:
                        msg = (("Failed to calculate reference ET values on %s. " % curdate) +
                               ("With input values:\n %s.\n" % str(rec)) +
                               ("Due to error: %s" % e))
                        raise PCSEError(msg)
        
                        # update record with ET values value convert to cm/day
                        rec[hdf5vars.index("e0")]  = E0/10.0;
                        rec[hdf5vars.index("es0")] = ES0/10.0;
                        rec[hdf5vars.index("et0")] = ET0/10.0;
                        
                        # If the record is incomplete, don't append!
                        if self._record_is_ok(rec):
                            ts.append(tuple(rec));
        except Exception as e:
            self.logger.error(str(e));
        finally:
            return ts;   
        
    def _record_is_ok(self, rec, roffset=0):
        for j in range(1, len(rec) - roffset):
            if rec[j] == None:
                return False;
        return True;
    
    def _check_Netcdf4rasters(self, timeref, envelope):
        try:
            # Compare the extent of the given envelope with that of the other one
            key1 = self.Netcdf4rasters.keys()[0]; # first key
            for key in self.Netcdf4rasters:
                nc4r = self.Netcdf4rasters[key];
                
                # Check the time reference
                times = nc4r.getVariables(Netcdf4Envelope2D.TIME);
                if (times.units != u'days since ' + str(timeref)):
                    fn = os.path.basename(nc4r.getFilePath());
                    raise PCSEError("Time in file " + fn + " not expressed as days since " + str(timeref));
                
                # Prepare to check the dates for which data are available in the various files
                start_day = date(self.available_years[0], 1, 1);
                end_day = date(self.available_years[-1], 12, 31);
                
                # Get the first and last day for which there are data in the file and check
                if (timeref + timedelta(times[0]) - datetime.combine(start_day, time(0,0,0))).days > 0: 
                    msg = "File " + fn + " does not have data for days as early as " + str(start_day.date());
                    raise PCSEError(msg);
                if (timeref + timedelta(times[-1]) - datetime.combine(end_day, time(0,0,0))).days < 0:
                    msg = "File " + fn + " does not have data for all days until " + str(end_day.date());
                    raise PCSEError(msg);
                
                # It was already checked that the file contains data wrt. a relevant variable. Now check the extent 
                # and sorting of coordinates
                if (key == key1): continue; 
                nvlp = nc4r.getEnvelope();            
                if not nvlp.hasSameExtent(envelope):
                    raise PCSEError("Netcdf4 files do not cover the same extent");
                if not nvlp.compareSorting(envelope):
                    raise PCSEError("Netcdf4 files do not have their coordinates sorted in the same way");
            return True;
        except Exception as e:
            self.logger.error(str(e));
            return False;
        
    def check_location(self, netcdf4raster, colIndex, rowIndex):
        try:
            # Get hold of all the data relevant for the given location, for the available years
            # Lookup lat-lon in the netCDF and check that latitudes and longitudes are the same!
            longitude, latitude = netcdf4raster.getEnvelope().getXandYfromIndices(colIndex, rowIndex); 
            latitudes = netcdf4raster.getVariables(Netcdf4Envelope2D.LAT);
            longitudes = netcdf4raster.getVariables(Netcdf4Envelope2D.LON);
            abs_diff = abs(latitudes[rowIndex] - latitude);
            assert abs_diff < 0.01, "Latitudes not equal: " + str(latitudes[rowIndex]) + " " + str(latitude);
            abs_diff = abs(longitudes[colIndex] - longitude)
            assert abs_diff < 0.01, "Longitudes not equal!: " + str(longitudes[rowIndex]) + " " + str(longitude);
            return True;
        except Exception as e:
            self.logger.error(str(e));
            return False;
        
    def _getAngstromCoeff(self, deglat):
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
        
    def _getFullPath(self, fname):
        key1 = self.Netcdf4rasters.keys()[0]; # first key
        nc4r = self.Netcdf4rasters[key1];
        path = os.path.dirname(nc4r.getFilePath());
        path = os.path.join(path, fname);
        return os.path.normpath(path);
    
class ncdf_weather(IsDescription):
    day = Time32Col(pos=1)
    tmax = Float32Col(pos=2)
    tmin = Float32Col(pos=3)
    temp = Float32Col(pos=4)
    rain = Float32Col(pos=5)
    irrad = Float32Col(pos=6)
    wind = Float32Col(pos=7)
    vap = Float32Col(pos=8)
    e0 = Float32Col(pos=9)
    es0 = Float32Col(pos=10)
    et0 = Float32Col(pos=11);
    
def getOrderedKeys(items):
    # Initialise
    result = [''] * len(items);
    
    # Result is zero-based but pos is not - correct
    maxpos = len(items) - 1;
    for key in items:
        value = items[key];
        maxpos = max(maxpos, value._v_pos);
    shift = maxpos - (len(items) - 1);
    
    # Now place the keys in the right place of the array 'result' 
    for key in items:
        value = items[key];  
        result[value._v_pos - shift] = key;
    return result;
