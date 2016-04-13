# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2015
from __future__ import division
from gridenvelope2d import GridEnvelope2D;
from basetiffraster import BaseTiffRaster
from libtiff import TIFF
from itertools import izip
from math import fabs
import numpy as np
import const
import os

class RowTiffRaster(GridEnvelope2D, BaseTiffRaster):
    "A raster represented by 2 files, with extensions 'tif' and 'tfw'"
    "This class can ONLY deal with tiff files which have 1 row per strip"
    "Different layers - e.g. RGB - are as planes: contiguously = chunky = interleaved"
    "or separately = per channel. More info: http://www.fileformat.info/format/tiff/egff.htm"
    "It means that the number of planes and the planar configuration determines the shape"
    "of the array written as bitmapped data, with dimensions image_depth, image_height, "
    "image_width and samples. E.g. in the case of rgb and contiguous configuration, the last"
    "dimension of the array is expected to be 3 and the field samples per pixel will also be 3"
    # Private attributes
    __mode = 'r'
    __datatype = const.FLOAT;
    currow = -1;
    __envelope = None;
    __bits_per_sample = 8
    __sample_format = 1
    __samples_per_pixel = 1
    __numpy_type = np.uint8
    __itemsize = 1
    __layer_size = 1
    
    # Predefine __ReadStrip with a dummy function
    def dummy(self, *args): pass;
    __ReadStrip = dummy(0, 0, 1)
    __WriteStrip = dummy(0, 0, 1)

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
                
    def __get_sample_format(self, arr):
        result = None
        
        # Not considered: SAMPLEFORMAT_VOID=4 and SAMPLEFORMAT_COMPLEXINT=5
        if arr.dtype in np.sctypes['float']:
            result = 3 #SAMPLEFORMAT_IEEEFP
        elif arr.dtype in np.sctypes['uint']+[np.bool]:
            result = 1 #SAMPLEFORMAT_UINT
        elif arr.dtype in np.sctypes['int']:
            result = 2 #SAMPLEFORMAT_INT
        elif arr.dtype in np.sctypes['complex']:
            result = 6 #SAMPLEFORMAT_COMPLEXIEEEFP
        else:
            raise NotImplementedError(`arr.dtype`)
        return result
    
    def set_numpy_type(self, type):
        self.__numpy_type = type
    
    def get_numpy_type(self):
        return self.datafile.get_numpy_type(self.__bits_per_sample, self.__sample_format)
    
    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=100, nodatavalue=-9999.0, byteorder='II', compression=1):
        # Initialise
        super(RowTiffRaster, self).open(mode);
        
        # If file does not exist and mode[0] = 'w', create it!
        if (mode[0] == 'w'):
            self.__mode = mode
            self.datafile = TIFF.open(self.folder + os.path.sep + self.name, mode='w');
            self.__envelope = GridEnvelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            self.datafile.SetField("ImageWidth", ncols)
            self.datafile.SetField("ImageLength", nrows)
            self.datafile.SetField("BitsPerSample", self.__bits_per_sample)
            self.datafile.SetField("SamplesPerPixel", self.__samples_per_pixel)
            self.datafile.SetField("RowsPerStrip", 1)
            self.datafile.SetField("PlanarConfig", 1) # contiguous
            self.datafile.SetField("Orientation", 1)  # top left
            self.datafile.SetField("PageNumber", 1, 1)
            self.datafile.SetField("FillOrder", 1) # MSB2LSB
            self.datafile.SetField("Compression", compression)
            super(RowTiffRaster, self).set_extra_tags()
            
            # Prepare to write per strip
            if compression == 1: # none
                self.__WriteStrip = self.datafile.WriteRawStrip
            else:
                self.__WriteStrip = self.datafile.WriteEncodedStrip
            self.writeheader()
            return True;
        else: 
            # Open the file as well as the header file
            if os.path.exists(os.path.join(self.folder, self.name)):            
                self.datafile = TIFF.open(os.path.join(self.folder, self.name), mode='r');
                self.readheader();
                
                # Check whether found values warrant further execution
                self.ncols = int(self.datafile.GetField("ImageWidth"))
                self.nrows = int(self.datafile.GetField("ImageLength"))
                if int(self.datafile.GetField("RowsPerStrip")) > 1:
                    msg = "Image file stores data with more than 1 row per strip. This class is unable to handle this"
                    raise ValueError(msg)
                
                # Process further information from the header file
                self.xll = self.xul;
                if self.ycoords_sort == 'DESC':
                    self.yll = self.yul - self.nrows * self.dy;
                else:
                    self.yll = self.yul + self.nrows * self.dy; 
                    
                # Prepare to read the file
                self.__bits_per_sample = self.datafile.GetField("BitsPerSample")
                self.__sample_format = self.datafile.GetField("SampleFormat")
                self.__samples_per_pixel = self.datafile.GetField("SamplesPerPixel")
                self.__numpy_type = self.datafile.get_numpy_type(self.__bits_per_sample, self.__sample_format)
                self.__itemsize = self.__bits_per_sample / 8
                self.__layer_size = self.ncols * self.__samples_per_pixel * self.__itemsize
                if self.__datatype == const.INTEGER:
                    self.nodatavalue = int(self.datafile.GetField("GDAL_NODATA"))
                else:
                    self.nodatavalue = float(self.datafile.GetField("GDAL_NODATA"))
                super(RowTiffRaster, self).get_extra_tags()

                if self.datafile.GetField("Compression") == 1: # none
                    self.__ReadStrip = self.datafile.ReadRawStrip
                else:
                    self.__ReadStrip = self.datafile.ReadEncodedStrip  
                return True;
            else: return False;
    
    def next(self, parseLine=True):
        # Is it possible to proceed? Otherwise generate StopIteration
        self.currow += 1;
        if (self.currow > self.nrows): raise StopIteration;
        if (self.currow > int(self.datafile.NumberOfStrips())): raise StopIteration;
        
        # Read the next row
        buf = np.zeros((self.ncols, self.__samples_per_pixel), self.__numpy_type)
        self.__ReadStrip(self.currow, buf.ctypes.data, int(self.__layer_size))
        return buf      
    
    def writenext(self, sequence_with_data):
        # Write the next data if possible, otherwise generate StopIteration
        # We cannot know whether exactly 1 row is included or not.
        # Is it possible to proceed? Otherwise generate StopIteration
        if self.currow == -1:
            sample_format = self.__get_sample_format(sequence_with_data)
            self.datafile.SetField("SampleFormat", sample_format)
        self.currow += 1;
        if (self.currow > self.nrows): raise StopIteration;
        if (self.currow > int(self.datafile.NumberOfStrips())): raise StopIteration;
        
        # Ok, try to write the data
        try:
            #data_sequence = np.ascontiguousarray(sequence_with_data)
            size = self.ncols * self.__samples_per_pixel * sequence_with_data.itemsize
            if len(sequence_with_data.shape) == 1 or sequence_with_data.shape[1] == 1:
                self.__WriteStrip(self.currow, sequence_with_data.ctypes.data, int(size))
            else:
                size = size * sequence_with_data.shape[0]
                sequence_with_data = np.ascontiguousarray(sequence_with_data)
                self.__WriteStrip(self.currow, sequence_with_data.ctypes.data, int(size))
            
            return True
        except StopIteration:
            raise StopIteration
        except ValueError as e:
            print str(e)
            raise ValueError
        except Exception as e:
            raise IOError(str(e)); 
          
    def readheader(self):
        # Header has 6 lines - without labels!
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
            if self.ycoords_sort == 'ASC':
                hf.write(str(self.dy) + "\n");
            else:
                hf.write(str(-1 * self.dy) + "\n");
            hf.write(str(self.xll + 0.5 * self.dx) + "\n");
            hf.write(str(self.yll + self.nrows * self.dy - 0.5 * self.dx) + "\n");
        except Exception, e:
            msg = "Header file " + hdrFilename + " could not be written in folder " + self.folder;
            raise IOError(msg + "(" + str(e) + ")");
        
    def reset(self):
        self.currow = -1;  
        
    def get_colormap(self):
        return self.datafile.GetField("ColorMap")
        
    def write_colormap(self, rgbTable=None, palette=None):
        # Palette is assumed to be a list with a length that is 3 times the no.
        # of colors, filled with r, g and b values alternating each other
        if palette != None:
            rgb_gen = izip(palette[0::3], palette[1::3], palette[2::3])
            cmap = list(rgb_gen)
            numcolors = 1 << self.__bits_per_sample
            rgbTable = np.zeros((3, numcolors), np.int)
            for k in range(len(cmap)):
                rgbTable[:, k] = cmap[k]
        else:
            if rgbTable == None:
                raise ValueError("No palette and no RGB table provided")
        self.datafile.SetField("PhotoMetric", 3)
        self.datafile.SetField("ColorMap", rgbTable)
    
    def close(self):
        if self.__mode[0] == 'w':
            self.datafile.WriteDirectory()
        super(RowTiffRaster, self).close()
        