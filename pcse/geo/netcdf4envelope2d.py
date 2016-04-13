from .gridenvelope2d import GridEnvelope2D;

class Netcdf4Envelope2D(GridEnvelope2D):
    # Constants
    LON = 'lon';
    LAT = 'lat';
    TIME = 'time'
    
    # TODO: here the names for the dimensions are hardcoded, but they can be retrieved
    # from the netCDF4 file 
    def __init__(self, ds):
        # Initialise
        LON = self.LON;
        LAT = self.LAT;
        
        # Make sure you can call the __init__ method of the base class
        if ds != None:
            # Assume it's open; get hold of the dimensions and variables
            self.__dataset = ds;
            _dims = ds.dimensions;
            _vars = ds.variables;
            x_range= self._readRange(_vars[LON]);
            y_range = self._readRange(_vars[LAT]);
            
            # Now get the necessary properties for initialising an envelope object
            dx = GridEnvelope2D._getStep(x_range[0], x_range[1], len(_dims[LON]));
            dy = GridEnvelope2D._getStep(y_range[0], y_range[1], len(_dims[LAT]));
            xll = min(_vars[LON]) - 0.5*dx;
            yll = min(_vars[LAT]) - 0.5*dy;
            GridEnvelope2D.__init__(self, len(_dims[LON]), len(_dims[LAT]), xll, yll, dx, dy);
        
            # In some netCDF files the longitudes and latitudes are stored differently
            if _vars[LON][0] > _vars[LON][-1]: self.xcoords_sort = 'DESC';
            if _vars[LAT][0] < _vars[LAT][-1]: self.ycoords_sort = 'ASC';
    
    def _readRange(self, varxy):
        # Argument is either x or y
        minxy = min(varxy);
        maxxy = max(varxy);
        return [minxy, maxxy];
    
    def getDimension(self, ds): 
        result = 2;
        if ds != None:
            result = len(ds.dimensions.keys());
        return result;
    