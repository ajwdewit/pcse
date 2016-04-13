import os.path
import re
import tables
from collections import Sequence
from tables.exceptions import NoSuchNodeError
from raster import Raster
from .gridenvelope2d import GridEnvelope2D

# Class for reading quick and dirty HDF5 format that can store weather data
# in a raster format in an efficient way, for fast access
class Hdf5Raster(GridEnvelope2D, Raster):
    "A raster represented by 2 files, with extensions 'h5' and 'hdr'"
    
    # Constants
    DATAFILEXT = "h5";
    HEADEREXT = ""
    
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
    datafile = None;
    currow = 0;
    
    def __init__(self, filepath):
        # Initialise super class instance
        Raster.__init__(self, filepath)
        # Retrieve the name from the filepath and assign - incl. extension
        self.name = os.path.basename(filepath);
        # Also derive the folder
        self.folder = os.path.dirname(filepath);
    
    def open(self, mode, ncols=1, nrows=1, xll=0., yll=0., cellsize=1., nodatavalue=-9999.0,
        dataset_name="dummy", group_prefix="row", table_prefix="col", index_format="04i", variables=[], units=[]):
        # Initialise
        super(Hdf5Raster, self).open(mode);
        fpath = os.path.join(self.folder, self.name);
        if (mode[0] == 'w'):
            # Open the file
            self.datafile = tables.open_file(fpath, 'w');
            
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
                # Retrieve the data attributes from the attributes in the file
                self.datafile = tables.open_file(fpath, 'r');
                self.readheader();
                GridEnvelope2D.__init__(self, self.ncols, self.nrows, self.xll, self.yll, self.cellsize, self.cellsize);
                
                return True;
            else: return False;   
    
    def readheader(self):
        # All information is stored in the file itself, i.e. as attributes
        # Read attributes and assign them. Assume that the file is open.
        f = self.datafile;
        try:
            self.ncols = f.root._v_attrs.ncols;
            self.nrows = f.root._v_attrs.nrows;
            self.xll = f.root._v_attrs.xllcorner;        
            self.yll = f.root._v_attrs.yllcorner;        
            self.cellsize = f.root._v_attrs.cellsize;        
            self.nodatavalue = f.root._v_attrs.nodata_value;
            self.root_contains = f.root._v_attrs.dataset_name; 
            self.group_prefix = f.root._v_attrs.group_prefix; 
            self.table_prefix = f.root._v_attrs.table_prefix;
            self.index_format = f.root._v_attrs.index_format;
            self.variables = ', '.join(f.get_node(f.root, "variables", classname='Array').read());
            self.units = ', '.join(f.get_node(f.root, "units", classname='Array').read());
        except Exception, e:
            raise IOError("Problem encountered while reading attributes (" + str(e) + ")");
        
    def getDataFileExt(self):
        result = self.DATAFILEXT;
        try:
            result = os.path.splitext(self.datafile)[1].strip('.');
        finally:
            return result;
        
    def writeheader(self):
        # All this information is stored in the hdf5 file itself - write them as attributes 
        try:
            f = self.datafile;  
            
            # Write attributes to file
            f.root._v_attrs.ncols = self.ncols;
            f.root._v_attrs.nrows = self.nrows;
            f.root._v_attrs.xllcorner     = self.xll;
            f.root._v_attrs.yllcorner     = self.yll;
            f.root._v_attrs.cellsize      = self.cellsize;
            f.root._v_attrs.nodata_value = self.nodatavalue;
            f.root._v_attrs.dataset_name = self.dataset_name;
            f.root._v_attrs.group_prefix = self.group_prefix;
            f.root._v_attrs.table_prefix = self.table_prefix;
            f.root._v_attrs.index_format = self.index_format;
            
            # Add arrays to the file - variables and units are already lists
            try:
                f.get_node(f.root, "variables", classname='Array')
                f.remove_node(f.root, "variables");
            except tables.NoSuchNodeError: 
                pass;
            f.create_array(f.root, "variables", self.variables);
            try:
                f.get_node(f.root, "units", classname='Array')
                f.remove_node(f.root, "variables");
            except tables.NoSuchNodeError: 
                pass;
            f.create_array(f.root, "units", self.units);
            f.flush();
        except Exception, e:
            msg = "Attributes could not be written to file: " + self.datafile.filename + "\n";
            raise IOError(msg + str(e));
        
    def next(self, parseLine=True):
        # Read the next data slice if possible, otherwise generate StopIteration
        result = None;
        try:
            if (self.currow > self.nrows): raise StopIteration;
            if parseLine:
                grp_name = (self.group_prefix + "%" + self.index_format) % self.currow; 
                grp = self.datafile.get_node(self.datafile.root, grp_name);
                result = [None] * self.ncols;
                for k in range(0, self.ncols):
                    tbl_name = (self.table_prefix + "%" + self.index_format) % k;
                    tbl = self.datafile.get_node(grp, tbl_name);
                    result[k] = tbl.read();   
            self.currow += 1;  # row index is zero-based!  
            return result;
        except:
            raise StopIteration; 
        
    def writenext(self, sequence_with_data, recordClass):
        # Write the next data if possible, otherwise generate StopIteration
        # We assume that exactly 1 row is included, with for each pixel of the
        # current row an array with records or at least a None value
        try:   
            # Check input
            assert isinstance(sequence_with_data, Sequence), "Given input not of the expected type!"
            assert len(sequence_with_data) == self.ncols, "Input array does not have the expected size!";
            msg = "Input class reference does not inherit from tables.IsDescription!"
            assert issubclass(recordClass, tables.IsDescription), msg
            
            # Assume that the group does not yet exist; note: row index is zero-based!
            if (self.currow >= self.nrows): raise StopIteration;
            filter1 = tables.filters.Filters(complevel=1, complib='blosc', fletcher32=True);
            grp_name = (self.group_prefix + "%" + self.index_format) % self.currow; 
            f = self.datafile;
            grp = f.create_group(f.root, grp_name, 'represents a row');
            
            # Now loop and add a table for each column
            for k in range(0, self.ncols):
                recs = sequence_with_data[k];
                if (recs != None) and (type(recs) is list):
                    # Create the table and aAdd the records to it
                    tbl_name = (self.table_prefix + "%" + self.index_format) % k;
                    tbl = f.create_table(grp, tbl_name, recordClass, expectedrows=len(recs), filters=filter1);
                    tbl.append(recs);
                    tbl.flush();
            self.currow += 1;
            return True;
        except Exception, e:
            print str(e);            
            raise StopIteration;
        
    def write(self, colIndex, recordList, recordClass):
        # TODO test this!
        if (recordList == None) or (not type(recordList) is list):
            raise ValueError("Records were not provided in the form of a list.");
        msg = "Input class reference does not inherit from tables.IsDescription!"
        assert issubclass(recordClass, tables.IsDescription), msg
        
        # Initialise - it is assumed that the instance has been moved to the intended row already
        f = self.datafile;
        grp_name = (self.group_prefix + "%" + self.index_format) % self.currow; 
        try: 
            grp = f.get_node(f.root, grp_name);    
        except NoSuchNodeError:
            grp = f.create_group(f.root, grp_name, 'represents a row');
            
        # We've reached the right group - now we have to somehow get hold of the right table
        tbl_name = (self.table_prefix + "%" + self.index_format) % colIndex;
        try:    
            # If the table already exists, delete it!
            tbl = f.del_node_attr(grp, tbl_name); 
            f.del_node_attr(grp, tbl_name);
        except NoSuchNodeError: pass;

        try:    
            # We can get hold of the right table - do so and add the records to it 
            filter1 = tables.filters.Filters(complevel=1, complib='blosc', fletcher32=True); 
            tbl = f.create_table(grp, tbl_name, recordClass, expectedrows=len(recordList), filters=filter1);  
            tbl.append(recordList);
            tbl.flush();
        except Exception as e:
            raise IOError(str(e));
    
    def flush(self):
        self.datafile.flush();    
        
    def close(self):
        if self.datafile:
            self.datafile.close();
            self.datafile = None;        
                
    def reset(self):
        super(Hdf5Raster, self).reset()
        