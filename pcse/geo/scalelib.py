# -*- coding: latin-1 -*-
"""ScaleLib - a Python library for scaling raster data"""
from geo.gridenvelope2d import GridEnvelope2D
from raster import Raster
from math import floor
try:
    import numpy as np
    import numpy.ma as ma
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

class ScaleManager():
    __initialised = False
    __rg = None
    ncols = 1
    nrows = 1
    __factor = 1
    __koffset = 0 # offset for column
    __ioffset = 0 # offset for row
    
    def __init__(self, rg, xll, yll, ncols, nrows, factor):
        # Check the inputs
        if not isinstance(rg, Raster):
            raise TypeError("Input is not a raster!")
        if not isinstance(rg, GridEnvelope2D):
            raise TypeError("Input is not a grid envelope!")
        if int(factor) != factor:
            raise ValueError("Factor should be an integer!")
        if int(nrows) != nrows:
            raise ValueError("Number of rows should be an integer!")
        if int(ncols) != ncols:
            raise ValueError("Number of columns should be an integer!")
        
        # Find out how many pixels the new LL is shifted relative to the LL of the input raster. The horizontal
        # and vertical shift should be a whole no. of pixels. Correct if the difference is not too much
        xoffset = (xll - rg.xll) / rg.cellsize
        yoffset = (yll - rg.yll) / rg.cellsize
        if (abs(round(xoffset) - xoffset) > 0.05) or (abs(round(yoffset) - yoffset) > 0.05):
            raise ValueError("Indicated LL corner does not coincide with grid point!")
        xll = rg.xll + round(xoffset) * rg.cellsize
        yll = rg.yll + round(yoffset) * rg.cellsize
        
        # Result raster should not exceed the boundaries of the input raster grid
        yul = rg.yll + rg.nrows * rg.cellsize
        A = (xll + factor * rg.cellsize * ncols > rg.xll + rg.ncols * rg.cellsize)
        B = (yll + factor * rg.cellsize * nrows > yul)
        C = (xll < rg.xll)     
        if A or B or C:
            raise ValueError("One of the corners of the intended result grid envelope is out of bounds!")
        
        # Simple assignments
        self.__rg = rg
        self.ncols = ncols
        self.nrows = nrows
        self.__factor = factor
        
        # If the output grid envelope does not exactly cover the input grid envelope
        # we may have to skip some of the first columns and rows 
        self.__koffset = int(round((xll - rg.xll) / rg.cellsize))
        self.__ioffset = int(round((yul - nrows * rg.cellsize * factor - yll) / rg.cellsize))
        for _ in range(self.__ioffset): rg.next()
        self.__initialised = True
        
    def __next_strip(self):
        # nrows and ncols are the length and width of the output raster
        result = None
        try:
            # 
            if HAS_NUMPY:
                if (self.__rg.dataformat == 'f'):
                    result = np.zeros((self.__factor, (self.__rg.ncols - self.__koffset)), dtype=np.float32)
                else: result = np.zeros((self.__factor, (self.__rg.ncols - self.__koffset)), dtype=np.int32) 
            else: 
                raise NotImplementedError("No implementation yet for Python without numpy")
            while True:
                row_in_strip = (self.__rg.currow - self.__ioffset) % self.__factor # also zero-based!
                curstrip = int(floor((self.__rg.currow - self.__ioffset) / self.__factor))
                if (curstrip >= int(self.nrows)): raise StopIteration
                row = self.__rg.next()
                result[row_in_strip, 0:(self.__rg.ncols - self.__koffset)] = row[self.__koffset:]
                if row_in_strip == self.__factor - 1:
                    return result 
        except Exception as e:
            print str(e)
    
    def next_aggregated(self, technique):
        # Initialise
        result = None
        try:
            if HAS_NUMPY:
                A = (self.__rg.dataformat == 'f')
                B = (technique == 'MEDIAN' and self.__factor % 2 == 0)
                C = (technique == 'MEAN')
                if A or B or C:
                    result = np.zeros((self.ncols), dtype=np.float32)
                else: result = np.zeros((self.ncols), dtype=np.int32) 
            else: 
                raise NotImplementedError("No implementation yet for Python without numpy")
             
            # Now aggregate
            curstrip = self.__next_strip()
            if curstrip != None:
                for k in range(0, self.ncols):
                    j = k * self.__factor
                    result[k] = self.__get_value(curstrip[0:self.__factor, j:j+self.__factor], technique)
            else: raise Exception("Unknown error")
            return result
        except StopIteration:
            raise StopIteration
        except Exception as e:
            print str(e)

    def __get_value(self, square, technique):
        # SUM — The sum (total) of the input cell values (default)
        # MAX — The largest value of the input cells
        # MEAN — The average value of the input cells 
        # MEDIAN — The median value of the input cells 
        # MIN — The smallest value of the input cells
        result = 0
        rel_diff = lambda a, b: (a - b) / b 
        mask = (abs(rel_diff(square, self.__rg.nodatavalue)) < 0.01)
        masq = ma.masked_array(square, mask)
        if not self.__has_data(masq):
            result = self.__rg.nodatavalue
        else:
            if technique == 'SUM':
                result = masq.sum()
            elif technique.startswith('MAX'):
                result = masq.max()
            elif technique.startswith('MIN'):
                result = masq.min()
            elif technique == 'MEAN':
                result = masq.mean()
            elif technique == 'MEDIAN':
                result = self.__get_median(masq.flatten())
        return result  
            
    def __has_data(self, sequence):
        result = True
        flatseq = sequence.flatten()
        if ma.count_masked(flatseq) == len(flatseq):
            result = False
        return result
    
    def __get_median(self, sequence):
        result = 0
        sequence.sort()
        if len(sequence) % 2 == 1:
            pos = int(floor(len(sequence) / 2))
            result = sequence[pos]
        else:
            pos = -len(sequence) / 2
            result = (sequence[pos-1] + sequence[pos]) / 2
        return result
        
    def close(self):
        self.__rg = None
        self.__initialised = False
