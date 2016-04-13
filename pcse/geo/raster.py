# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2015
import os;

# Abstract base class
class Raster(object):
    # Constants
    DATAFILEXT = "xxx";
    HEADEREXT = "xxw"
    datafile = None
    name = "dummy";
    folder = os.getcwd();
    nodatavalue = -9999.0;
    cellsize = 1;
    currow = 0
    
    def __init__(self, filepath):
        pass
    
    def open(self, mode, ncols=1, nrows=1, xll=0.0, yll=0.0, cellsize=1.0, nodatavalue=-9999.0):
        pass
    
    def readheader(self):
        pass

    def writeheader(self):
        pass
    
    def __iter__(self):
        return self; 
    
    def next(self, parseLine=True):
        pass

    def writenext(self, sequence_with_data):
        # input is sequence type - e.g. list, array.array or numpy.array
        pass
    
    @staticmethod
    def getDataFileExt(self):
        return self.DATAFILEXT;
    
    @staticmethod
    def getHeaderFileExt(self):
        return self.HEADEREXT 
    
    def close(self):
        if self.datafile:
            if hasattr(self.datafile, 'closed'):
                if not self.datafile.closed:
                    self.datafile.close()
            else:
                self.datafile.close()

    def reset(self):
        self.currow = 0; 
    
    