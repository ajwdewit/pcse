import os.path;
import re;
import tables
from .gridenvelope2d import GridEnvelope2D;

# Class for reading quick and dirty HDF5 format that can store weather data
# in a raster format in an efficient way, for fast access
class Hdf5Raster(GridEnvelope2D):
    "A raster represented by 2 files, with extensions 'h5' and 'hdr'"
    
    # Constants
    DATAFILEXT = "h5";
    HEADEREXT = "hdr";
    
    # Data attributes - assign some dummy values for the mean time
    name = "dummy.h5";
    folder = os.getcwd();
    cellsize = 1;
    nodatavalue = -9999.0;
    dataset_name = "dummy";
    group_prefix = "row";
    table_prefix = "col";
    index_format = "04i";
    variables = "temp, rain";
    units = "degrees Celsius, mm/day";
    
    # Private attributes
    __datafile = None;
    __currow = 0;
    
    def __init__(self, filepath):
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);
    
    def open(self, mode, ncols=1, nrows=1, xll=0, yll=0, cellsize=1, nodatavalue=-9999.0,
        dataset_name="dummy", group_prefix="row", table_prefix="col", index_format="04i", variables=[], units=[]):
        # Initialise
        fpath = os.path.join(self.folder, self.name);
        if (mode[0] == 'w'):
            # Open the file
            self.__datafile = tables.open_file(fpath, 'w');
            
            # Assign the data attributes 
            self.ncols = ncols;
            self.nrows = nrows;                    
            self.xll = xll;
            self.yll = yll;
            self.cellsize = cellsize;
            self.nodatavalue = nodatavalue;
            self.dataset_name = dataset_name;
            self.group_prefix = group_prefix;
            self.table_prefix = table_prefix;
            self.index_format = index_format;
            self.variables = variables;
            self.units = units;
            self.writeheader();
        else: 
            # If file does not exist, then ...
            if os.path.exists(fpath):
                # Retrieve the data attributes from the header file
                self.readheader();
                GridEnvelope2D.__init__(self, self.ncols, self.nrows, self.xll, self.yll, self.cellsize, self.cellsize);
                self.__datafile = tables.open_file(fpath, 'r');
                return True;
            else: return False;   
    
    def readheader(self):
        # Read header file and assign all attributes 
        pos = str.rfind(str(self.name), "." + self.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + self.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        fpath = os.path.join(self.folder, hdrFilename)
        if os.path.exists(fpath):
            hf = open(fpath, 'r');
            hl = hf.readline();
            self.ncols = int(hl.replace('ncols', '').strip());
            hl = hf.readline();
            self.nrows = int(hl.replace('nrows', '').strip());
            hl = hf.readline();
            self.xll = float(hl.replace('xllcorner', '').strip());        
            hl = hf.readline();
            self.yll = float(hl.replace('yllcorner', '').strip());        
            hl = hf.readline();
            self.cellsize = float(hl.replace('cellsize', '').strip());        
            hl = hf.readline();
            self.nodatavalue = float(hl.replace('NODATA_value', '').strip());
            hl = hf.readline();
            self.root_contains = hl.replace('dataset_name', '').strip();  
            hl = hf.readline();
            self.group_prefix = hl.replace('group_prefix', '').strip(); 
            hl = hf.readline();
            self.table_prefix = hl.replace('table_prefix', '').strip();
            hl = hf.readline();
            self.index_format = hl.replace('index_format', '').strip();
            hl = hf.readline();
            self.variables = re.sub(r'\s', '', hl.replace('variables', '')).split(',');
            hl = hf.readline();
            self.units = re.sub(r'\s', '', hl.replace('units', '')).split(',');
            hf.close();
        else: 
            msg = "Header file " + hdrFilename + " not found in folder " + self.folder;
            raise IOError(msg);
        
    @staticmethod
    def getDataFileExt(self):
        return self.DATAFILEXT;
    
    @staticmethod
    def getHeaderFileExt(self):
        return self.HEADEREXT; 
        
    def writeheader(self):
        # Write header file with all attributes 
        pos = str.rfind(str(self.name), "." + self.DATAFILEXT);
        if pos != -1: hdrFilename = self.name[0:pos] + "." + self.HEADEREXT
        else: raise ValueError("Invalid file name: " + self.name);
        try:
            # Open the file if it exists, otherwise create it
            fpath = os.path.join(self.folder + os.path.sep + hdrFilename);
            if os.path.exists(fpath):
                hf = open(fpath, 'w');
            else:
                hf = file(fpath, 'w');
   
            # Now write all the attributes
            hf.write("ncols         " + str(self.ncols) + "\n");
            hf.write("nrows         " + str(self.nrows) + "\n");
            hf.write("xllcorner     " + str(self.xll) + "\n");
            hf.write("yllcorner     " + str(self.yll) + "\n");
            hf.write("cellsize      " + str(self.cellsize) + "\n");
            hf.write("NODATA_value  " + str(self.nodatavalue) + "\n");
            hf.write("dataset_name  " + self.dataset_name + "\n");
            hf.write("group_prefix  " + self.group_prefix + "\n");
            hf.write("table_prefix  " + self.table_prefix + "\n");
            hf.write("index_format  " + self.index_format + "\n");
            hf.write("variables     " + ", ".join(self.variables) + "\n");
            hf.write("units         " + ", ".join(self.units) + "\n");
        except Exception, e:
            print e;
            msg = "Header file " + hdrFilename + " could not be written in folder " + self.folder;
            raise IOError(msg);
        
    def __iter__(self):
        return self;
        
    def next(self, parseLine=True):
        # Read the next data slice if possible, otherwise generate StopIteration
        result = None;
        try:
            if (self.__currow > self.nrows): raise StopIteration;
            if parseLine:
                grp_name = (self.group_prefix + "%" + self.index_format) % self.__currow; 
                grp = self.__datafile.get_node(self.__datafile.root, grp_name);
                result = [None] * self.ncols;
                for k in range(0, self.ncols):
                    tbl_name = (self.table_prefix + "%" + self.index_format) % k;
                    tbl = self.__datafile.get_node(grp, tbl_name);
                    result[k] = tbl.read();   
            self.__currow += 1;  # row index is zero-based!  
            return result;
        except:
            raise StopIteration; 
        
    def writenext(self, array_with_data, recordClass):
        # Write the next data if possible, otherwise generate StopIteration
        # We assume that exactly 1 row is included, with for each pixel of the
        # current row an array with records or at least a None value
        try:   
            # Check input
            assert type(array_with_data) is list, "Given input not of the expected type!"
            assert len(array_with_data) == self.ncols, "Input array does not have the expected size!";
            msg = "Input class reference does not inherit from tables.IsDescription!"
            assert issubclass(recordClass, tables.IsDescription), msg
            
            # Assume that the group does not yet exist; note: row index is zero-based!
            if (self.__currow >= self.nrows): raise StopIteration;
            filter1 = tables.filters.Filters(complevel=1, complib='blosc', fletcher32=True)
            grp_name = (self.group_prefix + "%" + self.index_format) % self.__currow; 
            f = self.__datafile;
            grp = f.create_group(f.root, grp_name, 'represents a row');
            
            # Now loop and add a table for each column
            rc_len = len(recordClass.columns);
            for k in range(0, self.ncols):
                recs = array_with_data[k];
                if (recs != None) and (type(recs) is list):
                    # Create the table and aAdd the records to it
                    tbl_name = (self.table_prefix + "%" + self.index_format) % k;
                    tbl = f.create_table(grp, tbl_name, recordClass, expectedrows=len(recs), filters=filter1);
                    tbl.append(recs);
                    tbl.flush();
            self.__currow += 1;
            return True;
        except Exception, e:
            print e;            
            raise StopIteration
        
    def flush(self):
        self.__datafile.flush();    
        
    def close(self):
        if self.__datafile:
            self.__datafile.close();
            self.__datafile = None;        
                
    def reset(self):
        self.__currow = 0; 