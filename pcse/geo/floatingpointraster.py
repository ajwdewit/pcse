import os.path;
import const;
import sys;
import types;
from .gridenvelope2d import GridEnvelope2D;

class FloatingPointRaster(GridEnvelope2D):
    "A raster represented by 2 files, with extensions 'flt' and 'hdr'"
    
    # Constants
    const.LSBFIRST = "LSBFIRST";
    const.DATAFILEXT = "flt";
    const.HEADEREXT = "hdr";
    const.BYTESPERCELL = 4;
    #const.epsilon = 0.00000001;
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy.flt";
    folder = os.getcwd();
    """
    ncols = 1;
    nrows = 1;
    xll = 0;
    yll = 0;
    """
    cellsize = 100;
    nodatavalue = -9999.0;
    byteorder = const.LSBFIRST;
    
    # Private attributes
    __datatype = const.FLOAT;
    __datafile = None;
    __currow = 0;
    
    def __init__(self, filepath, *datatype):
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);
        # Finally set the datatype
        if len(datatype) > 0:
            if (datatype[0] == const.INTEGER): 
                self.__datatype = const.INTEGER;
            else: 
                self.__datatype = const.FLOAT;        
        
    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=100, nodatavalue=-9999.0, byteorder=const.LSBFIRST):
        # If file does not exist and mode[0] = 'w', create it!
        if (mode[0] == 'w') and (not os.path.exists(self.folder + os.path.sep + self.name)):
            self.__datafile = file(self.folder + os.path.sep + self.name, 'w');
            GridEnvelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            return True;
        else:    
            # Open the file
            if os.path.exists(self.folder + os.path.sep + self.name):            
                self.__datafile = open(self.folder + os.path.sep + self.name, mode[0] + 'b'); 
                if (mode[0] == 'w'):
                    # Assign the data attributes 
                    self.ncols = ncols;
                    self.nrows = nrows;                    
                    self.xll = xll;
                    self.yll = yll;
                    self.cellsize = cellsize;
                    self.nodatavalue = nodatavalue;
                    self.byteorder = byteorder;
                    self.writeheader();
                else: 
                    # Retrieve the data attributes from the header file
                    self.readheader();
                GridEnvelope2D.__init__(self, self.ncols, self.nrows, self.xll, self.yll, self.cellsize, self.cellsize);
                return True;
            else: return False; 
          
    
    def readheader(self):
        # Read header file and assign all attributes 
        pos = str.rfind(self.name, "." + const.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + const.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        if os.path.exists(self.folder + os.path.sep + hdrFilename):
            hf = open(self.folder + os.path.sep + hdrFilename, 'r');
            hl = hf.readline();
            self.ncols = int(hl.replace('ncols', '').strip());
            hl = hf.readline();
            self.nrows = int(hl.replace('nrows', '').strip());
            hl = hf.readline();
            self.xll = float(hl.replace('xllcorner', '').strip());        
            hl = hf.readline();
            self.yll = float(hl.replace('yllcorner', '').strip());        
            hl = hf.readline();
            self.cellsize = float(hl.replace('cellsize', '').strip());        
            hl = hf.readline();
            if (self.__datatype == const.INTEGER): 
                self.nodatavalue = int(hl.replace('NODATA_value', '').strip());
            else: 
                self.nodatavalue = float(hl.replace('NODATA_value', '').strip());
            hl = hf.readline();
            self.byteorder = hl.replace('byteorder', '');         
            hf.close();
        else: 
            msg = "Header file " + hdrFilename + " not found in folder " + self.folder;
            raise IOError(msg);
    
    def __iter__(self):
        return self;
        
    def next(self):
        # Read the next row if possible, otherwise generate StopIteration
        try:
            self.__currow += 1;
            if (self.__currow > self.nrows): raise StopIteration;
            return self.__datafile.read(self.ncols * const.BYTESPERCELL);        
        except:
            raise StopIteration;       
    
    """
    def hasSameExtent(self, obj):
        if not isinstance(obj, FloatingPointRaster):
            return False;
        if self.ncols != obj.ncols:
            return False;
        if self.nrows != obj.nrows:
            return False;
        if abs(self.cellsize - obj.cellsize) > const.epsilon:
            return False;
        if abs(self.xll - obj.xll) > const.epsilon:
            return False;
        if abs(self.yll - obj.yll) > const.epsilon:
            return False;
        else:
            return True;
    """
        
    @staticmethod
    def getDataFileExt():
        return const.DATAFILEXT;
    
    @staticmethod
    def getHeaderFileExt():
        return const.HEADEREXT;    
    
    @staticmethod
    def getBytesPerCell():
        return const.BYTESPERCELL;
    
    def writeheader(self):
        # Write header file with all attributes 
        pos = str.rfind(self.name, "." + const.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + const.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        try:
            # Open the file if it exists, otherwise create it
            if os.path.exists(self.folder + os.path.sep + hdrFilename):
                hf = open(self.folder + os.path.sep + hdrFilename, 'w');
            else:
                hf = file(self.folder + os.path.sep + hdrFilename, 'w');
   
            # Now write all the attributes
            hf.write("ncols         " + str(self.ncols) + "\n");
            hf.write("nrows         " + str(self.nrows) + "\n");
            hf.write("xllcorner     " + str(self.xll) + "\n");
            hf.write("yllcorner     " + str(self.yll) + "\n");
            hf.write("cellsize      " + str(self.cellsize) + "\n");
            hf.write("NODATA_value  " + str(self.nodatavalue) + "\n");
            hf.write("byteorder     " + self.byteorder + "\n");
        except Exception, e:
            print e;
            msg = "Header file " + hdrFilename + " could not be written in folder " + self.folder;
            raise IOError(msg);
        
    def writenext(self, array_with_data):
        # Write the next data if possible, otherwise generate StopIteration
        # We cannot know whether exactly 1 row is included or not.
        try:          
            return self.__datafile.write(array_with_data);
        except Exception, e:
            print e;            
            raise StopIteration
        
    def close(self):
        if self.__datafile:
            if not self.__datafile.closed:
                self.__datafile.close();        
                
    def reset(self):
        self.__datafile.seek(0);
        self.__currow = 0;   
        