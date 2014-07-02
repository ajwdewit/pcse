from .netcdf4envelope2d import Netcdf4Envelope2D;
from .gridenvelope2d import GridEnvelope2D;
from netCDF4 import Dataset;
import os;

class Netcdf4Raster(Netcdf4Envelope2D):
    # Constants
    DATAFILEXT = 'nc4';
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy.nc4";
    folder = os.getcwd();
    __dataset = None;
    __varname = "";
    __currow = 0;
    cellsize = 1;
    nodatavalue = -9999.0;
    
    def __init__(self, filepath):
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);

    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=1, nodatavalue=-9999.0):
        # If file does not exist and mode[0] = 'w', create it!
        fpath = os.path.join(self.folder, self.name);
        if (mode[0] == 'w') and (not os.path.exists(fpath)):
            raise Exception("Writing of netCDF4 not implemented yet!");
            # Netcdf4Envelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            return False;
        else:
            if os.path.exists(fpath):  
                # Open the netCDF4 file          
                ds = Dataset(fpath, 'r');
                Netcdf4Envelope2D.__init__(self, ds);
                self.__dataset = ds;
                
                # Establish which variable is stored in it - assume it's only 1!
                diff = set(ds.variables) - set(ds.dimensions);
                if len(diff) > 0: self.__varname = diff.pop();
                return True;
            else: return False;    
            
    def readheader(self):
        pass;
    
    def __iter__(self):
        return self;
        
    def next(self, parseLine=True):
        # The netCDF4 format is not one that is read in a sequential manner. Is it a good idea
        # that this method is part of the interface?
        result = None;
        if self.__varname != "":
            if parseLine:
                result = self.getVariables(self.__varname)[:, self.__currow, :];
            self.__currow += 1;  # row index is zero-based! 
            return result;
    
    @staticmethod
    def getDataFileExt(self):
        return self.DATAFILEXT;
    
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
        self.__currow = 0; 
    
    def getEnvelope(self):
        # TODO: this cannot be solved by invoking super or so?
        return GridEnvelope2D(self.ncols, self.nrows, self.xll, self.yll, self.dx, self.dy);
        
    def getVariables(self, dimDescriptor):
        return self.__dataset.variables[dimDescriptor];
    
    def getFilePath(self):
        return self.__dataset.filepath();
    
    def getVariableName(self):
        return self.__varname;
        