import math;

# A proposed superclass for AsciiGrid
class GridEnvelope2D:
    # Constants
    epsilon = 0.0000001; 
    
    # Data attributes
    nrows = 1;
    ncols = 1;
    xll = 0.0;
    yll = 0.0;
    dx = 2.0;
    dy = 2.0;

    # We define a grid by means of a number of rows, columns, 
    # a lower left corner as well as by steps in x and y direction
    def __init__(self, nrows, ncols, xll, yll, dx, dy):
        self.nrows = nrows;
        self.ncols = ncols;
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

    def getRowAndColIndex(self, x, y):
        # Return the zero-based row and column numbers of the grid cell
        # which houses the given point; count rows from the top downwards
        eps = self.epsilon; # we'll use this constant to prevent negative indices
        k = int(round((x - self.getMinX() - 0.5*self.dx + eps) / self.dx));
        i = int(round((self.getMaxY() - 0.5*self.dy - y + eps) / self.dy));
        return i, k;
        
    def getNearestCenterPoint(self, x, y):
        # Rows are counted from the top downwards
        i, k = self.getRowAndColIndex(x, y);
        cx = self.xll + (k + 0.5) * self.dx;
        cy = self.getMaxY() - (i + 0.5) * self.dy;
        return (cx, cy);
    
    def hasSameExtent(self, obj):
        if (not isinstance(obj, GridEnvelope2D)): 
            return False;
        if self.ncols != obj.ncols:
            return False;
        if self.nrows != obj.nrows:
            return False;
        if abs(self.dy - obj.dy) > self.epsilon:
            return False;
        if abs(self.dx - obj.dx) > self.epsilon:
            return False;
        if abs(self.xll - obj.xll) > self.epsilon:
            return False;
        if abs(self.yll - obj.yll) > self.epsilon:
            return False;
        else:
            return True;        
    
    @staticmethod
    def getStep(mincoord, maxcoord, numcells):
        # Coordinates & number of cells either in X or in Y direction
        return math.ceil(maxcoord - mincoord) / numcells;

    