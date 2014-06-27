from .netcdf4envelope2d import Netcdf4Envelope2D;
from netCDF4 import Dataset;
import const;
import os;

class Netcdf4Raster(Netcdf4Envelope2D):
    # Constants
    const.DATAFILEXT = 'nc4';
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy.nc4";
    folder = os.getcwd();
    __dataset = None;
    
    def __init__(self, filepath):
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);

    def open(self, mode):
        # If file does not exist and mode[0] = 'w', create it!
        fpath = os.path.join((self.folder, os.path.sep + self.name));
        if (mode[0] == 'w') and (not os.path.exists(fpath)):
            raise Exception("Writing of netCDF4 not implemented yet!");
            # Netcdf4Envelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            return False;
        else:
            if os.path.exists(fpath):            
                self.__dataset = Dataset(fpath, 'r');
            
            else: return False;    
            
    def readheader(self):
        pass;
    
    def __iter__(self):
        return self;
        
    def next(self):
        # TODO: netCDF4 is not a file format that is read in a sequential manner
        # Is it appropriate that this method is part of the interface? Is it possible
        # to return more dimensional results - how is this facilitated?
        # Do we only offer to return rows or also columns if desired?
        pass;
    
    @staticmethod
    def getDataFileExt():
        return const.DATAFILEXT;
    
    def writeheader(self):
        # appropriate?
        raise Exception("Writing of netCDF4 not implemented yet!");
    
    def writenext(self, array_with_data):
        # Not implemented yet
       raise Exception("Writing of netCDF4 not implemented yet!");
        
    def close(self):
        if self.__dataset:
            self.__dataset.close(); 

    def reset(self):
        # TODO: is it appropriate that this method is part of the interface?
        pass; 