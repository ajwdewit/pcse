"""ClipLib - a Python library for clipping from rasters"""
from gridenvelope2d import GridEnvelope2D
from array import array
try:
    import numpy as np
    import numpy.ma as ma
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# Determines the minimum envelope around an area with the same raster values     
def get_envelope(raster, code):
    # Initialise
    result = None
    mink = raster.ncols - 1
    mini = raster.nrows - 1
    maxk = 0
    maxi = 0
    
    try:
        # Loop over the rows
        for i in range(raster.nrows):
            # Assume that the land mass of the country consists of 1 part
            line = raster.next()
            left = most_left(line, code)
            right = most_right(line, code)
            
            # maxi and maxk are the largest indices
            if left != -1: mink = min(left, mink)
            if right != -1: maxk = max(right, maxk)
            if (left != -1) and (mini == raster.nrows - 1): mini = i
            if (left == -1) and (mini != raster.nrows - 1) and (maxi == 0): maxi = i - 1
                    
        # Now calculate the dimension of the rectangle
        if maxi == 0: maxi = raster.nrows - 1
        nrows = maxi - mini + 1
        ncols = maxk - mink + 1
        xll = raster.xll + mink * raster.dx
        yll = raster.yll + (raster.nrows - maxi - 1) * raster.dy
        result = GridEnvelope2D(ncols, nrows, xll, yll, raster.dx, raster.dy)
    except Exception as e:
        print str(e)
    finally:
        raster.reset()   
        return result

def is_clippable(grid, code):
    result = False
    try:
        # Now loop over the rows
        for _ in range(grid.nrows):
            line = grid.next()
            indices = np.where(line == code)[0]
            if len(indices) > 0: 
                result = True
                break
    finally:
        grid.reset()
        return result    

def rel_diff(a, b):
    return (a - b) / b

# Method get_mask to mask out all pixels which are not part of the area
# Input is overall zone grid, envelope and area code. Optional is the grid for which
#  a mask should be made. If included, pixels with missing values are added to the mask 
def get_mask(zonegrid, nvlp, value, valuegrid=None):
    # Initialise
    shape = (nvlp.nrows, nvlp.ncols)
    result = np.zeros(shape, dtype=np.int)
    mink = int(round((nvlp.xll - zonegrid.xll) / nvlp.dx))
    maxk = mink + nvlp.ncols - 1
    maxi = zonegrid.nrows - int(round((nvlp.yll - zonegrid.yll) / nvlp.dy)) - 1
    mini = maxi - nvlp.nrows + 1
    
    try:
        # Now loop over the rows
        for i in range(zonegrid.nrows):
            if i < mini: 
                zonegrid.next(False)
                if valuegrid != None: valuegrid.next(False)
                continue
            if i > maxi: break
            zline = zonegrid.next()

            # Use numpy functionality to get the mask:
            if isinstance(zline, (np.ndarray, np.generic)):
                ztest = zline[mink:maxk+1]
            else:
                ztest = np.array(zline[mink:maxk+1])
            ztest.shape = (nvlp.ncols)
                
            # If the value grid was specified, then we're interested in missing values
            vtest = None
            if valuegrid != None: 
                # Use numpy functionality here too
                vline = valuegrid.next()
                if isinstance(vline, (np.ndarray, np.generic)):
                    vtest = vline[mink:maxk+1]
                else:
                    vtest = np.array(vline[mink:maxk+1])
              
            # Now test  
            if vtest == None:
                result[i - mini, :] = (ztest != value) 
            else:
                blntmp = (abs(rel_diff(vtest, valuegrid.nodatavalue)) < 0.01)
                blntmp.shape = ztest.shape
                result[i - mini, :] = np.logical_or((ztest != value), blntmp) 
    except Exception as e:
        print str(e) 
    finally:
        zonegrid.reset()
        if valuegrid != None: valuegrid.reset()
        return result
        
def get_datatype(dataformat):
    result = np.bool
    if dataformat in ('h', 'H', 'i', 'I'):
        result = np.int
    elif dataformat == 'l':
        result = np.longlong
    elif dataformat == 'f':
        result = np.float32
    elif dataformat =='d':
        result = np.float64 
    return result

# Method get_masked_clip_result
def get_masked_clip_result(datagrid, nvlp, mask, initvalue=0):
    # The datagrid is a raster type
    # Extract data from the datagrid and return it as a 2-dimensional array
    shape = (nvlp.nrows, nvlp.ncols)
    datatype = get_datatype(datagrid.dataformat)
    if initvalue == 0: result = np.zeros(shape, dtype=datatype)
    
    # Now calculate which part of the datagrid we need
    # Assume that the datagrid coincides exactly with the "country pixels" 
    mink = int(round((nvlp.xll - datagrid.xll) / nvlp.dx))
    maxk = mink + nvlp.ncols - 1
    maxi = datagrid.nrows - int(round((nvlp.yll - datagrid.yll) / nvlp.dy)) - 1
    mini = maxi - nvlp.nrows + 1
    
    # Loop over the rows - assume that the file position is zero
    for i in range(datagrid.nrows):
        if i < mini: 
            datagrid.next(False)
            continue
        if i > maxi: break
        line = datagrid.next()
        arrline = line[mink:maxk+1]
        if isinstance(arrline, array):
            arrline = np.array(arrline) 
        arrline.shape = (nvlp.ncols)
        result[i - mini, :] = arrline
        
    # Wind up
    datagrid.reset()
    result = ma.masked_array(result, mask)
    return result
        
def most_left(line, value):
    result = -1
    if not HAS_NUMPY:
        for k in range(len(line)):
            if line[k] == value:
                result = k
                break
    else:
        test = np.where(line==value)[0]
        if len(test) > 0: result = test[0]
        else: result = -1
    return result   

def most_right(line, value):
    result = -1
    if not HAS_NUMPY:
        for k in reversed(range(len(line))):
            if line[k] == value:
                result = k
                break
    else:
        test = np.where(line==value)[0]
        if len(test) > 0: result = test[-1]
        else: result = -1
    return result 

def get_clipped_data(country_rg, data_rg, countrycode):
    # Initialise
    result = (None, None)
    try:
        # Check that the country code is a value in the country grid
        if not is_clippable(country_rg, countrycode):
            raise Warning("Given country code not found!")
        
        # Get the envelope
        nvlp = get_envelope(country_rg, countrycode) 
        if nvlp == None:
            raise Exception("Not able to construct an envelope!")
        
        # Get a mask that coincides with the envelope and the array with the relevant values 
        mask = get_mask(country_rg, nvlp, countrycode, data_rg)
        result = (nvlp, get_masked_clip_result(data_rg, nvlp, mask)) # 2nd element is an array
    except Exception as e:
        print str(e)
    finally:
        return result

def join_envelopes(nvlps):
    # Initialise
    xll = 180
    yll = 90
    xur = -180
    yur = -90
    dx = 0.1
    dy = 0.1
    eps = 0.0001;
    
    # Check whether the envelopes can be joined
    for i, nvlp in zip(range(len(nvlps)), nvlps):
        if i == 0:
            dx = nvlp.dx
            dy = nvlp.dy
        else:
            if abs(nvlp.dx - dx) > eps or abs(nvlp.dy - dy) > eps:
                msg = "Cell sizes of the grids are not the same!"
                raise  Exception(msg)
    
    for nvlp in nvlps:  
        xll = min(xll, nvlp.xll)
        yll = min(yll, nvlp.yll)
        xur = max(xur, nvlp.getMaxX())
        yur = max(yur, nvlp.getMaxY())
        
    ncols = int(round((xur - xll) / dx))
    nrows = int(round((yur - yll) / dy))
    result = GridEnvelope2D(ncols, nrows, xll, yll, dx, dy)
    return result
           
# Method create_raster that is to return an envelope and an array with data
# Input is a minimum envelope and a masked array with values
def create_raster(countrydata, envelope=None, nodatavalue=None, debug="n"):
    # Assume that the countrydata are tuples with envelope and masked arrays
    if envelope == None:
        total_nvlp = join_envelopes([cd[0] for cd in countrydata])
    else:
        total_nvlp = envelope
        
    # Assume that all arrays hold the same data type
    shape = (total_nvlp.nrows, total_nvlp.ncols)
    datagrid = countrydata[0][1] 
    mask = np.ones(shape, dtype=np.int)
    total_array = np.ma.zeros(shape, datagrid.dtype)
    if nodatavalue != None: total_array.fill(nodatavalue)
    total_array.mask = mask
    
    # Now loop over the countries
    j = 0
    for cd in countrydata:
        # Figure out where we'll have to write
        (nvlp, clip_result) = cd
        xoffset = int(round((nvlp.xll - total_nvlp.xll) / total_nvlp.dx))
        yoffset = int(round((total_nvlp.getMaxY() - nvlp.getMaxY()) / total_nvlp.dy))
        
        # Now assign the clip_result - assume it's also a numpy array
        if len(debug) > 0 and debug[0].lower() == "y":
                j += 1
                print "About to add values for country number " + str(j)
        total_array = assign_slice2d(total_array, clip_result, (yoffset,yoffset+nvlp.nrows), (xoffset,xoffset+nvlp.ncols))
    
    total_array.mask = ma.nomask
    result = (total_nvlp, total_array)
    return result
    
def assign_slice2d(arr_dest, arr_src, seqidx0, seqidx1):
    # Check whether both input arrays are numpy masked arrays and the last arguments tuples of length 2
    if not isinstance(arr_dest, ma.masked_array) or not isinstance(arr_src, ma.masked_array):
        raise ValueError("At least one of the input arrays is not a numpy masked array as expected")
    if not isinstance(seqidx0, tuple) or not isinstance(seqidx1, tuple):
        raise ValueError("Last 2 arguments are not tuples, as expected!")
    if len(seqidx0) != 2 or len(seqidx1) != 2:
        raise ValueError("Last 2 arguments are not both of length 2, as expected!")
    
    # Check whether the area indicated by the sequence indices is the same as the dimensions of arr_src
    len0 = seqidx0[1] - seqidx0[0] # rows
    len1 = seqidx1[1] - seqidx1[0] # cols
    if arr_src.shape != (len0, len1):
        raise ValueError("Size of given source array does not correspond with given sequence indices!")
    
    # Now construct an index for the destination array
    dest_shp = arr_dest.shape
    dest_idx = np.ma.array(range(dest_shp[0] * dest_shp[1]), dtype=np.int)
    dest_idx.shape = dest_shp
    
    # Prepare an index for the destination array that is masked for those places that are masked in the source array
    dest_idx = dest_idx[seqidx0[0]:seqidx0[1], seqidx1[0]:seqidx1[1]]
    dest_idx.mask = arr_src.mask
    dest_idx = dest_idx.flatten().compressed()
    result = arr_dest.flatten()
    
    # Prepare an index for the source array that is masked for those places that are masked in the source array
    src_shp = arr_src.shape
    src_idx = np.ma.array(range(src_shp[0] * src_shp[1]), dtype=np.int)
    src_idx.shape = arr_src.shape 
    src_idx.mask = arr_src.mask
    src_idx = src_idx.flatten().compressed() 
    
    # Now use the prepared indices to do the assignment
    arr_src = arr_src.flatten()
    result[dest_idx] = arr_src[src_idx]
    result.shape = dest_shp
    return result
