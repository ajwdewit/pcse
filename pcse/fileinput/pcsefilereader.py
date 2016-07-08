# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import os, sys
import inspect
import textwrap

class PCSEFileReader(dict):
    """Reader for parameter files in the PCSE format.
    
    This class is a replacement for the `CABOFileReader`. The latter can be
    used for reading parameter files in the CABO format, however this format
    has rather severe limitations: it only supports string, integer, float
    and array parameters. There is no support for specifying parameters with
    dates for example (other then specifying them as a string).
    
    The `PCSEFileReader` is a much more versatile tool for creating parameter
    files because it leverages the power of the python interpreter for
    processing parameter files through the `execfile` functionality in python.
    This means that anything that can be done in a python script can also be
    done in a PCSE parameter file.

    :param fname: parameter file to read and parse
    :returns: dictionary object with parameter key/value pairs.

    *Example*
    
    Below is an example of a parameter file 'parfile.pcse'. Parameters can
    be defined the 'CABO'-way, but also advanced functionality can be used by
    importing modules, defining parameters as dates or numpy arrays and even
    applying function on arrays (in this case `np.sin`)::

        \"\"\"This is the header of my parameter file.
        
        This file is derived from the following sources
        * dummy file for demonstrating the PCSEFileReader
        * contains examples how to leverage dates, arrays and functions, etc.
        \"\"\"
        
        import numpy as np
        import datetime as dt
        
        TSUM1 = 1100
        TSUM2 = 900
        DTSMTB = [ 0., 0.,
                   5., 5.,
                  20., 25.,
                  30., 25.]
        AMAXTB = np.sin(np.arange(12))
        cropname = "alfalfa"
        CROP_START_DATE = dt.date(2010,5,14)

    Can be read with the following statements::
    
        >>>fileparameters = PCSEFileReader('parfile.pcse')
        >>>print fileparameters['TSUM1']
        1100
        >>>print fileparameters['CROP_START_DATE']
        2010-05-14
        >>>print fileparameters
        PCSE parameter file contents loaded from:
        D:\UserData\pcse_examples\parfile.pw
        
        This is the header of my parameter file.

        This file is derived from the following sources
        * dummy file for demonstrating the PCSEFileReader
        * contains examples how to leverage dates, arrays and functions, etc.
        DTSMTB: [0.0, 0.0, 5.0, 5.0, 20.0, 25.0, 30.0, 25.0] (<type 'list'>)
        CROP_START_DATE: 2010-05-14 (<type 'datetime.date'>)
        TSUM2: 900 (<type 'int'>)
        cropname: alfalfa (<type 'str'>)
        AMAXTB: [ 0.          0.84147098  0.90929743  0.14112001 -0.7568025
          -0.95892427  -0.2794155   0.6569866   0.98935825  0.41211849
          -0.54402111 -0.99999021] (<type 'numpy.ndarray'>)
        TSUM1: 1100 (<type 'int'>)
    """
    
    def __init__(self, fname):
        dict.__init__(self)
        
        # Construct full path to parameter file and check file existence
        cwd = os.getcwd()
        self.fname_fp = os.path.normpath(os.path.join(cwd, fname))
        if not os.path.exists(self.fname_fp):
            msg = "Could not find parameter file '%s'" % self.fname_fp
            raise RuntimeError(msg)

        # compile and execute the contents of the file
        bytecode = compile(open(self.fname_fp).read(), self.fname_fp, 'exec')
        exec(bytecode, {}, self)
        
        # Remove any members in self that are python modules
        keys = list(self.keys())
        for k in keys:
            if inspect.ismodule(self[k]):
                self.pop(k)
        
        # If the file has a header (e.g. __doc__) store it.
        if "__doc__" in self:
            header = self.pop("__doc__")
            if len(header) > 0:
                self.header = header
                if self.header[-1] != "\n":
                    self.header += "\n"
        else:
            self.header = None

    def __str__(self):
        printstr = "PCSE parameter file contents loaded from:\n"
        printstr += "%s\n\n" % self.fname_fp
        if self.header is not None:
            printstr += self.header
        for k in self:
             r = "%s: %s (%s)" % (k, self[k], type(self[k]))
             printstr += (textwrap.fill(r, subsequent_indent="  ") + "\n")
        return printstr
