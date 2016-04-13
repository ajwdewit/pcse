# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2015
from gridenvelope2d import GridEnvelope2D;
from raster import Raster
from math import fabs;
import const, os;
from tifffile import TiffFile;
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False 

class GeotiffRaster(GridEnvelope2D, Raster):
    "A raster represented by 2 files, with extensions 'tif' and 'tfw'"
    
    # Constants
    DATAFILEXT = "tif";
    HEADEREXT = "tfw"; # WORLD_EXT?
    __BYTESPERCELL = 4;
    #const.epsilon = 0.00000001;
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy.tif";
    folder = os.getcwd();
    nodatavalue = -9999.0;
    byteorder = 'II'; # Little endian
    roty = 0.0;
    rotx = 0.0;
    
    # Private attributes
    __datatype = const.FLOAT;
    #datafile = None;
    currow = -1;
    __envelope = None;
    __image = None;
    
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
    
    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=100, nodatavalue=-9999.0, byteorder='II'):
        # Initialise
        super(GeotiffRaster, self).open(mode);
        
        # If file does not exist and mode[0] = 'w', create it!
        if (mode[0] == 'w') and (not os.path.exists(self.folder + os.path.sep + self.name)):
            self.datafile = file(self.folder + os.path.sep + self.name, 'w');
            self.__envelope = GridEnvelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            return True;
        else: 
            # Open the file
            if os.path.exists(self.folder + os.path.sep + self.name):            
                self.datafile = TiffFile(self.folder + os.path.sep + self.name, multifile=False);
                if (mode[0] == 'w'):
                    # Assign the data attributes. TODO distinguish between dx and dy
                    dx = cellsize;
                    dy = cellsize; 
                    self.ncols = ncols;
                    self.nrows = nrows;                    
                    self.xul = xll;
                    if self.ycoords_sort == 'DESC':
                        self.yul = yll + nrows * dy;
                    else:
                        self.yul = yll - nrows * dy;
                    self.dx = dx;
                    self.dy = dy;
                    self.writeheader();
                else: 
                    # Retrieve some data attributes from the header file
                    self.readheader();
                    self.xll = self.xul;
                        
                    # Retrieve other ones from the TIFF file itself 
                    if self.datafile.byteorder == '<': self.byteorder = 'II';
                    else: self.byteorder = 'MM';
                    self.__image = self.datafile.asarray();
                    if len(self.__image.shape) > 3:
                        raise ValueError("Not sure how to handle data with more than 3 dimensions!");
                    axes = self.datafile.series[0]['axes'];
                    posx = axes.index('X')
                    posy = axes.index('Y')
                    posc = 3 - posx - posy
                    self.ncols = self.__image.shape[posx]; 
                    self.nrows = self.__image.shape[posy];
                    if HAS_NUMPY:
                        self.__image = np.rollaxis(self.__image, posc, 1)
                    else:
                        raise ImportError("Numpy not found on this system")
                    msg = "Data not properly dimensioned"
                    if not self.__image.shape[1] == self.nrows: raise Exception(msg)
                    if not self.__image.shape[2] == self.ncols: raise Exception(msg)
                    if self.ycoords_sort == 'DESC':
                        self.yll = self.yul - self.nrows * self.dy;
                    else:
                        self.yll = self.yul + self.nrows * self.dy; 
                    try:
                        page = self.datafile.pages[0];
                        if 'gdal_nodata' in page.tags: self.nodatavalue = float(page.gdal_nodata);
                    except: pass; 
                self.__envelope = GridEnvelope2D.__init__(self, self.ncols, self.nrows, self.xll, self.yll, self.dx, self.dy);
                return True;
            else: return False;
    
    def readheader(self):
        # header has 6 lines - without labels!
        sign = lambda x: (1, -1)[x<0];
        pos = str.rfind(str(self.name), "." + self.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + self.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        if os.path.exists(self.folder + os.path.sep + hdrFilename):
            # Adapt the following so that it accounts also for rotated mapsheets
            hf = open(self.folder + os.path.sep + hdrFilename, 'r');
            hl = hf.readline();
            self.dx = float(hl.strip());
            hl = hf.readline();
            self.roty = float(hl.strip());
            hl = hf.readline();
            self.rotx = float(hl.strip());
            eps = 0.0001;
            if abs(self.rotx)>eps or abs(self.roty)>eps:
                raise NotImplementedError("Cannot handle rotated mapsheets yet!")
            hl = hf.readline();            
            self.dy = fabs(float(hl.strip()));
            if sign(float(hl.strip())) == 1.0: self.ycoords_sort = 'ASC';
            hl = hf.readline();
            self.xul = float(hl.strip()) - 0.5 * self.dx;
            hl = hf.readline();
            self.yul = float(hl.strip()) + 0.5 * self.dy;
            hf.close();
        else: 
            msg = "Header file " + hdrFilename + " not found in folder " + self.folder;
            raise IOError(msg);
        
    def next(self, parseLine=True):
        # Read the next row if possible, otherwise generate StopIteration
        result = None;
        try:
            self.currow += 1;
            if (self.currow >= self.nrows): raise StopIteration;
            if parseLine:
                if len(self.__image.shape) == 1:
                    result = self.__image[self.currow];  
                else:
                    result = self.__image[:, self.currow]
            return result     
        except:
            raise StopIteration; 
    
    def writeheader(self):
        # Write header file with all attributes 
        pos = str.rfind(str(self.name), "." + self.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + self.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        try:
            # Open the file if it exists, otherwise create it
            if os.path.exists(self.folder + os.path.sep + hdrFilename):
                hf = open(self.folder + os.path.sep + hdrFilename, 'w');
            else:
                hf = file(self.folder + os.path.sep + hdrFilename, 'w');
   
            # Now write all the attributes
            hf.write(str(self.dx) + "\n");
            eps = 0.0001;
            if abs(self.rotx)>eps or abs(self.roty)>eps:
                raise NotImplementedError("Cannot handle rotated mapsheets yet!")
            hf.write(str(self.rotx) + "\n");
            hf.write(str(self.roty) + "\n");
            if self.ycoords_sort == 'DESC':
                hf.write(str(self.dy) + "\n");
            else:
                hf.write(str(-1 * self.dy) + "\n");
            hf.write(str(self.xul + 0.5 * dx) + "\n");
            hf.write(str(self.yul + 0.5 * dy) + "\n");
        except Exception, e:
            msg = "Header file " + hdrFilename + " could not be written in folder " + self.folder;
            raise IOError(msg + "(" + str(e) + ")");
    
    def writenext(self, sequence_with_data):
        # Write the next data if possible, otherwise generate StopIteration
        # We cannot know whether exactly 1 row is included or not.
        pass;
                
    def reset(self):
        self.currow = -1;   
    
    def get_colormap(self):
        try:
            page = self.datafile.pages[0];
            return page.color_map
        except:
            print "Unable to get color map" 
            
    def get_numpy_type(self):
        if HAS_NUMPY:
            result = np.int8
            if self.__image != None:
                result = self.__image.dtype
            return result
        else:
            ImportError("Numpy not found on this system")
 