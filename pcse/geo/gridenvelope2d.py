import math;
import const;

# Superclass for AsciiGrid and other classes
class GridEnvelope2D:
    # Constants
    #epsilon = 0.0000001;
    const.FLOAT = 'f';
    const.INTEGER = 'i'; 
    const.XLLCORNER = "XLLCORNER";    
    const.YLLCORNER = "YLLCORNER";
    const.epsilon = 0.00000001;
    
    # Data attributes
    nrows = 1;
    ncols = 1;
    xll = 0.0;
    yll = 0.0;
    dx = 2.0;
    dy = 2.0;
    xcoords_sort = 'ASC';  # ascending
    ycoords_sort = 'DESC'; # descending

    # We define a grid by means of a number of rows, columns, 
    # a lower left corner as well as by steps in x and y direction
    def __init__(self, ncols, nrows, xll, yll, dx, dy):
        self.ncols = ncols;
        self.nrows = nrows;
        self.xll = xll;
        self.yll = yll;
        self.dx = dx;
        self.dy = dy;

    def getDimension(self): 
        return 2;

    # The following 4 functions are about the extent of the envelope    
    def getMinX(self):
        # X-coordinate of the lower left corner
        return self.xll;
    
    def getMinY(self):
        # Y-coordinate of the lower left corner
        return self.yll;
    
    def getMaxX(self):
        # X-coordinate of the upper right corner
        return self.xll + self.ncols * self.dx;
    
    def getMaxY(self):
        # Y-coordinate of the upper right corner
        return self.yll + self.nrows * self.dy;

    def getColAndRowIndex(self, x, y):
        # Return the zero-based row and column numbers of the grid cell
        # which houses the given point; count rows either from the bottom upwards (ASC)
        # or from the top downwards (DESC), depending on the way the coordinates are sorted
        eps = const.epsilon; # we'll use this constant to prevent negative indices
        if self.xcoords_sort == 'ASC':
            k = int(round((x - self.xll - 0.5*self.dx + eps) / self.dx));
        else:
            k = int(round((self.getMaxX() - x - 0.5*self.dy + eps) / self.dy)); # TODO: check
        if self.ycoords_sort == 'ASC':
            i = int(round((y - self.yll - 0.5*self.dy + eps) / self.dy));
        else:
            i = int(round((self.getMaxY() - y - 0.5*self.dy + eps) / self.dy)); # TODO: check
        return k, i;
    
    def getXandYfromIndices(self, k, i):
        # TODO implement for xcoords_sort == 'DESC' and ycoords_sort == 'ASC'
        lon = float('inf');
        lat = float('inf');
        if self.xcoords_sort == 'ASC':
            lon = self.getMinX() + k*self.dx + 0.5*self.dx;
        if self.ycoords_sort == 'DESC':
            lat = self.getMaxY() - i*self.dy - 0.5*self.dy;
        return lon, lat;
        
    def getNearestCenterPoint(self, x, y):
        # Rows are counted from the top downwards
        k, i = self.getColAndRowIndex(x, y);
        if self.xcoords_sort == 'ASC':
            cx = self.xll + (k + 0.5) * self.dx;
        else:
            cx = self.getMaxX() - (k + 0.5) * self.dx;
        if self.ycoords_sort == 'ASC':
            cy = self.yll + (i + 0.5) * self.dy;
        else:
            cy = self.getMaxY() - (i + 0.5) * self.dy;
        return (cx, cy);
    
    def hasSameExtent(self, obj):
        if (not isinstance(obj, GridEnvelope2D)): 
            return False;
        if self.ncols != obj.ncols:
            return False;
        if self.nrows != obj.nrows:
            return False;
        if abs(self.dy - obj.dy) > const.epsilon:
            return False;
        if abs(self.dx - obj.dx) > const.epsilon:
            return False;
        if abs(self.xll - obj.xll) > const.epsilon:
            return False;
        if abs(self.yll - obj.yll) > const.epsilon:
            return False;
        else:
            return True;        
    
    def compareSorting(self, obj):
        result = True;
        if (not isinstance(obj, GridEnvelope2D)): 
            result = False;
        if self.xcoords_sort != obj.xcoords_sort:
            result = False;
        if self.ycoords_sort != obj.ycoords_sort:
            result = False;
        return result; 
    
    @staticmethod
    def _getStep(mincoord, maxcoord, numcells):
        # Coordinates & number of cells either in X or in Y direction
        return math.ceil(maxcoord - mincoord) / numcells;
    
    def isWithinExtent(self, x, y):
        A = (x >= self.getMinX());
        B = (x <= self.getMaxX());
        C = (y >= self.getMinY());
        D = (y <= self.getMaxY());
        return (A and B and C and D);


    