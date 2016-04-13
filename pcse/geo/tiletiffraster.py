# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), March 2016
from __future__ import division
from gridenvelope2d import GridEnvelope2D;
from basetiffraster import BaseTiffRaster
from libtiff import TIFF, libtiff
import numpy as np
import os
import const
from math import floor, ceil, sqrt, fabs

class TileTiffRaster(GridEnvelope2D, BaseTiffRaster):
    "A raster represented by 2 files, with extensions 'tif' and 'tfw'"
    "This class can deal with tiff files of which the data are stored in tiles"
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
    __image = None; # contains current strip
    __bits_per_sample = 8
    __sample_format = 1
    __samples_per_pixel = 1
    __numpy_type = np.uint8
    __itemsize = 1
    __layer_size = 1
    __tile_width = 1
    __tile_length = 1
    __ntiles = 1
    __nstrips = 1 
    
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
    
    def set_numpy_type(self, atype):
        self.__numpy_type = atype
    
    def get_numpy_type(self):
        return self.datafile.get_numpy_type(self.__bits_per_sample, self.__sample_format)
    
    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=100, nodatavalue=-9999.0, byteorder='II', compression=1):
        # Initialise
        super(TileTiffRaster, self).open(mode);
        
        # If file does not exist and mode[0] = 'w', create it!
        if (mode[0] == 'w'):
            self.__mode = mode
            self.datafile = TIFF.open(self.folder + os.path.sep + self.name, mode='w');
            self.__envelope = GridEnvelope2D.__init__(self, ncols, nrows, xll, yll, cellsize, cellsize);
            self.datafile.SetField("ImageWidth", ncols)
            self.datafile.SetField("ImageLength", nrows)
            self.datafile.SetField("BitsPerSample", self.__bits_per_sample)
            self.datafile.SetField("SamplesPerPixel", self.__samples_per_pixel)
            super(TileTiffRaster, self).set_extra_tags()
             
            # Data are organised into square tiles. Let each tile be about 8K bytes
            bits_per_pixel = self.__samples_per_pixel * self.__bits_per_sample
            pixels_per_tile = max(int(floor(8 * 8000 / bits_per_pixel)), 1)
            self.__tile_width = floor(sqrt(pixels_per_tile))
            self.__tile_length = self.__tile_width
            self.__ntiles = int(ceil(ncols / self.__tile_width) * ceil(nrows / self.__tile_length))
            self.__nstrips = int(self.ntiles / ceil(self.ncols / self.__tile_width))
            self.datafile.SetField("TileWidth", self.__tile_width)
            self.datafile.SetField("TileLength", self.__tile_length)
            self.datafile.SetField("PlanarConfig", 1) # contiguous
            self.datafile.SetField("Orientation", 1)  # top left
            self.datafile.SetField("PageNumber", 1, 1)
            self.datafile.SetField("FillOrder", 1) # MSB2LSB
            self.datafile.SetField("Compression", compression)
            self.writeheader()
            shape = (self.__tile_length * self.ncols, self.__samples_per_pixel)
            self.__image = np.zeros(shape, self.__numpy_type)
            return True;
        else: 
            # Open the file as well as the header file
            if os.path.exists(self.folder + os.path.sep + self.name):            
                self.datafile = TIFF.open(self.folder + os.path.sep + self.name, mode='r');
                self.readheader();
                
                # Check whether found values warrant further execution
                self.ncols = int(self.datafile.GetField("ImageWidth"))
                self.nrows = int(self.datafile.GetField("ImageLength"))
                self.__tile_width = int(self.datafile.GetField("TileWidth"))
                self.__tile_length = int(self.datafile.GetField("TileLength"))
                self.__ntiles = libtiff.TIFFNumberOfTiles(self.datafile).value
                
                # Tiles can be joined to form strips: ceil(ncols / tile_width) in number, with height = tile_length
                # Those strips can be joined to form the image: ceil(nrows / tile_length) in number 
                msg = "Number of tiles not in accordance with tile and image dimensions!"
                self.__nstrips = int(ceil(self.nrows / self.__tile_length))
                num_tiles_per_strip = int(ceil(self.ncols / self.__tile_width))
                assert self.__ntiles == self.__nstrips * num_tiles_per_strip, msg
                planar_config = self.datafile.GetField("PlanarConfig")
                if (planar_config > 1):
                    raise NotImplementedError("Not yet able to deal with data organised in separate planes")
                if self.__datatype == const.INTEGER:
                    self.nodatavalue = int(self.datafile.GetField("GDAL_NODATA"))
                else:
                    self.nodatavalue = float(self.datafile.GetField("GDAL_NODATA"))
                super(TileTiffRaster, self).get_extra_tags()
                
                # Process further information from the header file
                self.xll = self.xul;
                if self.ycoords_sort == 'DESC':
                    self.yll = self.yul - self.nrows * self.dy;
                else:
                    self.yll = self.yul + self.nrows * self.dy; 
                    
                # Prepare to read the file (strip by strip and under the hood tile by tile)
                self.__bits_per_sample = self.datafile.GetField("BitsPerSample")
                self.__sample_format = self.datafile.GetField("SampleFormat")
                self.__samples_per_pixel = self.datafile.GetField("SamplesPerPixel")
                self.__numpy_type = self.datafile.get_numpy_type(self.__bits_per_sample, self.__sample_format)
                self.__itemsize = self.__bits_per_sample / 8
                shape = (self.__tile_length * self.ncols, self.__samples_per_pixel)
                self.__image = np.zeros(shape, self.__numpy_type)
                return True;
            else: return False;
    
    def get_tag(self, name):
        return self.datafile.get_tag_name(name)
    
    def next(self, parseLine=True):
        # Is it possible to proceed? Otherwise generate StopIteration
        result = None;
        self.currow += 1;
        try:
            if (self.currow >= self.nrows): raise StopIteration;
        
            # Read a new strip when necessary
            row_in_strip = self.currow % self.__tile_length # also zero-based!
            curstrip = int(floor(self.currow / self.__tile_length)) 
            if curstrip >= self.__nstrips: raise StopIteration;
            if row_in_strip == 0:
                # Are we dealing with one plane or with more? What configuratin?
                # self.datafile.GetField("PlanarConfig", 1)) 
                if curstrip == self.__nstrips-1:
                    # Last strip
                    length = self.nrows - self.__tile_length * curstrip
                    self.__layer_size = (self.ncols) * length * (self.__samples_per_pixel) * (self.__itemsize) 
                    self.__image = np.zeros((length, self.ncols, self.__samples_per_pixel), dtype=self.__numpy_type) 
                    #.resize((last_length, self.ncols, self.__samples_per_pixel)) 
                else:
                    length = self.__tile_length
                    self.__layer_size = (self.ncols) * length * (self.__samples_per_pixel) * (self.__itemsize) 
                    

                # Before trying to read, reset the buffer           
                self.__image.fill(0.0)
                self.__ReadStrip(curstrip, self.__image, int(self.__layer_size))
                self.__image = self.__image.reshape(length, self.ncols, self.__samples_per_pixel)
          
            # Read the next row
            result = self.__image[row_in_strip, :]
            return result     
        except StopIteration:
            raise StopIteration;   
        except Exception as e:
            print str(e)
    
    def writenext(self, sequence_with_data):
        raise NotImplementedError("Not implemented yet")
    
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
    
    def __ReadStrip(self, strip, buf, size):
        result = False
        try:
            num_tiles_per_strip = int(ceil(self.ncols / self.__tile_width))
            numpy_type = self.datafile.get_numpy_type(self.__bits_per_sample, self.__sample_format)
            if (strip == self.__nstrips-1): 
                length = self.nrows - (strip*self.__tile_length)
            else:
                length = self.__tile_length
            buf = buf.reshape(length, self.ncols, self.__samples_per_pixel)
            for k in range(num_tiles_per_strip): 
                if (k == num_tiles_per_strip-1):
                    # We only need part of the tile because we are on the edge
                    width = self.ncols - (num_tiles_per_strip-1)*self.__tile_width
                else:
                    width = self.__tile_width
                tmp_buf = np.ascontiguousarray(np.zeros((self.__tile_length, self.__tile_width), numpy_type))
                seq = libtiff.TIFFReadTile(self.datafile, tmp_buf.ctypes.data, k*self.__tile_width, strip*self.__tile_length, 0, 0)
                if seq != None:
                    start = k*self.__tile_width
                    buf[0:length, start:start+width, 0] = tmp_buf[0:length, 0:width]
            result = True
        except Exception as e:    
            print str(e)
        finally:
            return result
    
    def __WriteStrip(self,strip, buf, size):
        raise NotImplementedError("Not implemented yet")
        
    def get_colormap(self):
        return self.datafile.GetField("ColorMap")
        
    def write_colormap(self, rgbTable=None, palette=None):
        pass
    
    def close(self):
        if self.__mode[0] == 'w':
            self.datafile.WriteDirectory()
        super(TileTiffRaster, self).close()
    
    