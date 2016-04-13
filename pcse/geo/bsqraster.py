# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), June 2015
import os
import stat
import struct
from math import fabs
from array import array
from bandraster import BandRaster
try:
    import numpy
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

class BsqRaster(BandRaster):
    "A raster represented by 2 files, with extensions 'bsq' and 'bsq.hdr'"
    
    # TODO: change some function calls so that it's more obvious how to use "with" statements
    # Use of such statements are essential when files are written to disk!
    
    # Data attributes - assign some dummy values for the mean time
    DATAFILEXT = 'bsq'
    HEADEREXT = "hdr"
    name = "dummy.bsq"
    folder = os.getcwd();
    nbands = 3 #default
    nbits = 8  #default
    dataformat = "h" # Default data format (2-byte signed short int)
    WORLDEXT = "bqw"
    
    def __init__(self, filepath, *dataformat):
        # Initialise super class instance
        BandRaster.__init__(self, filepath, dataformat[0])
        self.byteorder = self.INTEL
        self.pixeltype = self.UNSIGNEDINT
        self.currow = -1; 

    def open(self, mode, ncols=1, nrows=1, nbands=3, xll=0, yll=0, cellsize=100, nodatavalue=256):
        result = super(BsqRaster, self).open(mode, ncols, nrows, nbands, xll, yll, cellsize, nodatavalue)
        if (mode[0] == 'w'):
            return result
        else:
            # Open the file
            if os.path.exists(os.path.join(self.folder, self.name)):            
                # Check given format string is valid
                bytesperpix = 2 #default
                try:
                    bytesperpix = struct.calcsize(self.dataformat)
                except:
                    raise ValueError, "Supplied data format " + str(self.dataformat) + " is invalid"
                # end try
                
                # Check file size matches with size attributes
                fileinfo = os.stat(os.path.join(self.folder, self.name))
                filesize = fileinfo[stat.ST_SIZE]
                if (filesize == 0) and (mode[0] == 'w'):
                    print "Empty BSQ file found. I'm going to overwrite it ..."
                else:
                    checknum = (((filesize / float(self.nbands)) / float(self.nrows)) / float(bytesperpix)) / self.ncols
                    if checknum != 1:
                        if fabs(checknum - 1) < 0.00003:
                            raise ValueError, "File size and size calculated from attributes only match approximately"
                        else:
                            raise ValueError, "File size and supplied attributes do not match at all!"
            
                # Open the file for reading in binary mode
                try:
                    self.datafile = open(os.path.join(self.folder, self.name), mode[0] + "b")
                except:
                    msg = "Failed to open BSQ file " + os.path.join(self.folder, self.name)
                    raise IOError(msg)          
                return True;
            else: return False; 
    
    def next(self, parseLine=True):
        super(BsqRaster, self).next()

        # Read the next row if possible, otherwise generate StopIteration
        try:
            self.currow += 1;
            if (self.currow > self.nrows): raise StopIteration
            
            # Get the size in bytes for the given data format
            line = []
            if not parseLine:
                return line
            
            # Idea is to return the i-th line with values from each band. We're dealing with a file
            # that has its i-th line scattered over blocks of data each representing only 1 band
            itemsize = struct.calcsize(self.dataformat)
            rowsize = self.ncols * itemsize 
            blocksize =  self.nrows * rowsize
            bandstartpos = []
            for j in range(self.nbands):
                 bandstartpos.append(j * blocksize)   
            
            if not HAS_NUMPY:
                # The following is for the case without numpy: each line should be
                # a list of length self.ncols, with arrays of length self.nbands 
                for _ in range(0, self.ncols):
                    line.append(array(self.dataformat, self.nbands*[0]))
            else:
                line = numpy.zeros((self.ncols, self.nbands))
            
            # Read the right line from each block - even though this may be slow
            for bandnum in range(0, self.nbands):
                self.datafile.seek(bandstartpos[bandnum] + self.currow * rowsize)
                buffer = self.datafile.read(rowsize)
                dataitems = struct.unpack_from(self.ncols*self.dataformat, buffer)
                if not HAS_NUMPY:
                    for pixnum in range(self.ncols):
                        line[pixnum][bandnum] = dataitems[pixnum]
                    # end for
                else:
                    line[:, bandnum] = dataitems
            # end for
                
            return line            
                   
        except StopIteration:
            raise StopIteration
        except Exception, e:
            raise Exception(e)
    
    def writenext(self, sequence_with_data):
        super(BsqRaster, self).writenext(sequence_with_data)
        # TODO: test!
        # Write the next data if possible, otherwise generate StopIteration
        # Assume that the input sequence is a sequence of length self.ncols
        # with 1-dimensional values
        try:            
            # Perform a number of checks
            if not self._is_sequence(sequence_with_data):
                raise ValueError("Input value is not a sequence!")
            if len(sequence_with_data) != self.ncols:
                raise ValueError("Input sequence has not got the expected length")

            # We're writing 1 line for the current band. Assign the data to the right data structure
            itemsize = struct.calcsize(self.dataformat)
            buffer = bytearray(self.ncols * itemsize)
            line = sequence_with_data[0:self.ncols]
            struct.pack_into(self.ncols*self.dataformat, buffer, 0, *line)
            self.datafile.write(buffer)
            self.datafile.flush()
            
            return True
        except StopIteration:
            raise StopIteration
        except ValueError as e:
            print str(e)
            raise ValueError
        except Exception as e:
            raise IOError(str(e))
    
    def reset(self):
        self.currow = -1
        self.datafile.seek(0)
    
        