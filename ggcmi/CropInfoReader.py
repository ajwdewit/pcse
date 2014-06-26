import os;
from pcse.pcse.geo import GridEnvelope2D, AsciiGrid;
from netCDF4 import Dataset;
from pcse.pcse.fileinput import CABOFileReader;
from pcse.pcse.netcdf4 import NetcdfEnvelope2D;
from pcse.pcse.exceptions import PCSEError;
from datetime import datetime, date;
from numpy import ma;

class CropInfoProvider():
    regions = {"Europe": GridEnvelope2D(70, 136, -11., 36., 0.5, 0.5),
               "World": GridEnvelope2D(360, 720, -180., -90., 0.5, 0.5),
               "Tropics": GridEnvelope2D(94, 720, -180., -23.5, 0.5, 0.5)};
    
    # CABO name - filename
    # TODO add more elements for different continents ???
    _crop_info_sources = [("Barley", "World", "BAR301.CAB", 1),
            ("Cassava", "World", "CASSAVA.W41", 2),
            ("Groundnuts", "World", "GR_NUT.W41", 1),
            ("Maize", "Europe", "MAG203.CAB", 1), 
            ("Millet", "World", "MILLET.W41", 1),
            ("Potatoes", "World",  "POT702.CAB", 2),
            ("Pulses", "World", "PIGEOPEA.W41", 1),
            ("Rapeseed", "World", "RAP1002.CAB", 1),
            ("Rice", "World", "RIC501.CAB", 1),
            ("Rye", "World", "RYE.W41", 1),
            ("Sorghum", "World", "SORGHUM.W41", 1),
            ("Soybeans", "World", "SOYBEAN.W41", 1),
            ("Sunflower", "World", "SUN1101.CAB", 1),
            ("Wheat", "World", "WWH105.CAB", 1)
            ];
            # ("Maize", "Tropics", "MAIZE.W41", 1),
            
    # File - path is relative to files in subdirectory OtherInputs\GrowingSeason_vs1_24
    landmask_grid = r"../../geodata/glob_landmask_resampled.flt";
    
    # Netcdf dataset with planting and harvest dates
    _ds = None;
    _cropdata = None;    
    _envelope = None;
    
    def __init__(self, crop, watersupply="rf", fpath=None):
        # match the given crop and region with a tuple from the list
        _fn = self._crop_info_sources[0][2];
        for cropname, region, filename, group_no in self._crop_info_sources:
            # for the mean time, don't worry about the region
            if cropname == crop:
                _fn = filename;
                break;
        
        # read parameters from the filename stored in _crop_info_sources
        _fn = os.path.join(fpath, "OtherInputs", 'CROPD', _fn);
        self._cropdata = CABOFileReader(_fn);
    
        # Prepare to read planting and harvest data from the netcdf file
        _fn = crop + "_" + watersupply + "_growing_season_dates_v1.24.nc4";
        _fn = os.path.join(fpath, "OtherInputs", "GrowingSeason_vs1_24", _fn);

        try:
            # Open the file
            self._ds = Dataset(_fn, 'r');
            self._envelope = NetcdfEnvelope2D.getEnvelope(self._ds);
        except Exception as e:
            fn = os.path.basename(self._ds.filepath());
            raise PCSEError("An error occurred while opening file " + fn + " (" + str(e) + ")");
    
    @staticmethod
    def getCropGroup(aCropName):
        result = 0;
        for cropname, region, filename, group_no in CropInfoProvider._crop_info_sources:
            if cropname == aCropName:
                result = group_no;
                break;
        return result;
    
    @staticmethod
    def getCrops():
        result = {};
        i = 1;
        for cropname, region, filename, group_no in CropInfoProvider._crop_info_sources:
            result[i] = cropname;
            i = i + 1;
        return result;
    
    def getCropData(self):
        return self._cropdata;   
    
    def getSeasonDates(self, longitude, latitude):
        start_day = -99;
        end_day = -99;
        try:
            # Check that the netCDF file is loaded
            if self._ds == None: raise PCSEError("file is no more open.");
            
            # Check that the given latitude and longitude is within extent of the netCDF4 file
            if not self._envelope.isWithinExtent(longitude, latitude):
                raise PCSEError("Given lat-lon coordinates are beyond the borders of the file.");
            k, i = self._envelope.getColAndRowIndex(longitude, latitude);
            
            # Check that the found indices are really linked to the given lat-lon
            LON = NetcdfEnvelope2D.LON;
            LAT = NetcdfEnvelope2D.LAT;
            msg = " coordinate for this index not as expected: "
            assert self._ds.variables[LAT][i] == latitude, "Y" + msg + str(i);
            assert self._ds.variables[LON][k] == longitude, "X" + msg + str(k);
            
            # Look for the variables indicating the season
            planting_day = 'planting day';
            harvest_day  = 'harvest day';
            if (planting_day in self._ds.variables) and (harvest_day in self._ds.variables):
                # Take into account that the variables are made out of masked arrays
                arr_elem = self._ds.variables[planting_day][i, k];
                if str(arr_elem) != '--': start_day = int(arr_elem);
                arr_elem = self._ds.variables[harvest_day][i, k];
                if str(arr_elem) != '--': end_day = int(arr_elem);
            return start_day, end_day;
        except Exception as e:
            fn = os.path.basename(self._ds.filepath());
            raise PCSEError("An error occurred while reading file " + fn + " (" + str(e) + ")");
            
    def getTimerData(self, start_day, end_day, year): 
        result = None;
        try:
            # Prepare the timer data
            result = {};
            result['CAMPAIGNYEAR'] = year;
            if end_day > start_day:
                # Assume that start and end time are within the same year
                crop_end_date = datetime.strptime(str(year) + ' ' + str(end_day), '%Y %j').date();
                result['START_DATE'] = date(year, 1, 1);
                result['END_DATE'] = date(year, 12, 31);
            else:
                # Not within 1 year; try to start the season 90 days before the crop start date
                crop_end_date = datetime.strptime(str(year + 1) + ' ' + str(end_day), '%Y %j').date();
                tmpday = max(end_day, start_day - 90);
                result['START_DATE'] = datetime.strptime(str(year) + ' ' + str(tmpday), '%Y %j').date();
                result['END_DATE'] = datetime.strptime(str(year + 1) + ' ' + str(tmpday-1), '%Y %j').date();
            crop_start_date = datetime.strptime(str(year) + ' ' + str(start_day), '%Y %j').date();
            result['CROP_START_DATE'] = crop_start_date;        
            result['CROP_END_DATE'] = crop_end_date     
            result['CROP_START_TYPE'] = 'sowing';
            result['CROP_END_TYPE'] = 'harvest'; 
            result['MAX_DURATION'] = 300;
            return result;
        except Exception as e:
            raise PCSEError("An error occurred while preparing the timer data: " + str(e));
        
    def close(self):
        if self._ds != None: 
            self._ds.close();
            self._ds = None;
            
    def getExtent(self):
        return self._envelope;
        
    def _getFullPath(self, fname):
        path = os.path.dirname(self._ds.filepath());
        result = os.path.join(path, fname);
        result = os.path.normpath(result);
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
        if (k == r.ncols): k = k - 1;
        return line[int(k)];   
    
    def get_landmask(self, longitude, latitude):
        fpath = self._getFullPath(self.landmask_grid);
        value = self._get_value_from_grid(longitude, latitude, fpath);
        return (value == 1);
    
    
            