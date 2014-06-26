import const;
import os.path;
import array;
from .gridenvelope2d import GridEnvelope2D;

class AsciiGrid(GridEnvelope2D):
    "A raster represented by an ASCII file, with extension 'asc'"
    
    # Constants
    const.FILEXT = "asc";     
    const.MAXDIGITSPERCELL = 7;
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy." + const.FILEXT;
    folder = os.getcwd();
    """
    ncols = 1;
    nrows = 1;
    xll = 0.0;
    yll = 0.0;
    """
    cellsize = 100.0;
    nodatavalue = -9999.0;
    datatype = const.FLOAT;    
    
    # Private attributes
    __datafile = None;
    __currow = 0;
    __mode = 'r';
    __digitspercell = 7;
    
    def __init__(self, filepath, *datatype):
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);
        # Finally set the datatype
        if len(datatype) > 0:
            if (datatype[0] == const.INTEGER): 
                self.datatype = const.INTEGER;
            else: 
                self.datatype = const.FLOAT;
                
    def open(self, mode, ncols=1, nrows=1, xll=0.0, yll=0.0, cellsize=100.0, nodatavalue=-9999.0):       
        # If file does not exist and mode[0] = 'w', create it!
        if (mode[0] == 'w') and (not os.path.exists(self.folder + os.path.sep + self.name)):
            self.__datafile = file(self.folder + os.path.sep + self.name, 'w');
            self.__mode = mode;
            GridEnvelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            return True;
        else:    
            # Open the file
            if os.path.exists(self.folder + os.path.sep + self.name):            
                self.__datafile = open(self.folder + os.path.sep + self.name, mode[0]); 
                if (mode[0] == 'w'):
                    # Assign the data attributes 
                    self.ncols = ncols;
                    self.nrows = nrows;                    
                    self.xll = xll;
                    self.yll = yll;
                    self.cellsize = cellsize;
                    self.nodatavalue = nodatavalue;
                    self.writeheader();
                else: 
                    # File is open - retrieve the data attributes from the header of the file
                    self.readheader();
                    
                    # Also find out how many digits per cell were used - assume it's constant
                    pos = self.__datafile.tell();
                    line = self.__datafile.readline();
                    self.__digitspercell = ((1 + len(line)) / self.ncols) - 1;
                    self.__datafile.seek(pos);  # return to first line with data
                GridEnvelope2D.__init__(self, self.ncols, self.nrows, self.xll, self.yll, self.cellsize, self.cellsize);
                return True;
            else: return False;
            
    def readheader(self):
        # Assume that the file is open; read header of the file and assign all attributes 
        if (self.__datafile != None):
            if (not self.__datafile.closed):            
                hl = self.__datafile.readline();
                self.ncols = int(hl.replace('ncols', ''));
                hl = self.__datafile.readline();
                self.nrows = int(hl.replace('nrows', ''));
                hl = self.__datafile.readline();
                self.xll = float(hl.replace('xllcorner', ''));        
                hl = self.__datafile.readline();
                self.yll = float(hl.replace('yllcorner', ''));        
                hl = self.__datafile.readline();
                self.cellsize = float(hl.replace('cellsize', ''));        
                hl = self.__datafile.readline();
                if (self.datatype == const.INTEGER): 
                    self.nodatavalue = int(hl.replace('NODATA_value', ''));
                else: 
                    self.nodatavalue = float(hl.replace('NODATA_value', ''));
            else: 
                msg = "File " + self.name + " not found in folder " + self.folder;
                raise IOError(msg);                
    
    def __iter__(self):
        return self;
    
    def next(self, splitline=True):
        # Read the next row if possible, otherwise generate StopIteration
        # Assume that the header lines have been read and are correct wrt. ncols and nrows
        result = None;
        try:
            if (self.__datafile != None):
                if (not self.__datafile.closed):
                    self.__currow += 1;
                    if (self.__currow > self.nrows): 
                        raise StopIteration("Attempt to move beyond last row.");
                    
                    # Allocate a new array with ncols of the right type
                    if (self.datatype == const.INTEGER):
                        result = array.array('l', self.ncols * [self.nodatavalue]);
                    else:
                        result = array.array('f', self.ncols * [self.nodatavalue]);

                    # Now fill the array - first translate whitespace into space
                    rawline = self.__datafile.readline();
                    if splitline:
                        i = 0;                    
                        for x in rawline.split(): 
                            if (i < self.ncols):
                                if (self.datatype == const.INTEGER):
                                    result[i] = int(x);
                                else:
                                    result[i] = float(x);
                            i = i + 1;
                    return result;
                else: raise StopIteration("Attempt to read raster data from a closed file.");
            else: raise StopIteration("Attempt to read raster data from an unassigned file.")
        except Exception, e:
            print "Error: " + str(e);
            raise StopIteration;
    
    """       
    def hasSameExtent(self, obj):
        if not isinstance(obj, AsciiGrid):
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
    def getFileExt():
        return const.FILEXT;
        
    def writeheader(self):
        # Assume that the file is open; write header of the file with all attributes 
        if (self.__datafile != None):
            if (not self.__datafile.closed):
                try:
                    self.__datafile.write("ncols         " + str(self.ncols).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                    self.__datafile.write("nrows         " + str(self.nrows).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                    self.__datafile.write("xllcorner     " + str(self.xll).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                    self.__datafile.write("yllcorner     " + str(self.yll).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                    self.__datafile.write("cellsize      " + str(self.cellsize).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                    self.__datafile.write("NODATA_value  " + str(self.nodatavalue).rjust(const.MAXDIGITSPERCELL + 1) + "\n");
                except Exception, e:
                    print e;
                    msg = "Header lines could not be written to file " + self.name + " in folder " + self.folder;
                    raise IOError(msg);
        
    def writenext(self, array_with_data):
        # Write the next line if possible, otherwise generate StopIteration
        # We assume that exactly 1 row is included.
        try:          
            for i in range(0, self.ncols):
                s = str(array_with_data[i]).rjust(const.MAXDIGITSPERCELL + 1);
                self.__datafile.write(s);                
            return self.__datafile.write("\n");
        except Exception, e:
            print e;            
            raise StopIteration
    
    def flush(self):
        self.__datafile.flush();
    
    def close(self):
        if self.__datafile:
            if not self.__datafile.closed:
                self.__datafile.close();  
                
    def reset(self):
        self.__datafile.seek(0);
        if (self.__mode[0] == 'r'):
            self.readheader();
        self.__currow = 0;   
        