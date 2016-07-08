# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import re

from ..exceptions import PCSEError

class XYPairsError(PCSEError):
    pass

class LengthError(PCSEError):
    pass

class DuplicateError(PCSEError):
    pass
    
class CABOFileReader(dict):
    """Reads CABO files with model parameter definitions.
    
    The parameter definitions of Wageningen crop models are generally
    written in the CABO format. This class reads the contents, parses
    the parameter names/values and returns them as a dictionary.

    :param fname: parameter file to read and parse
    :returns: dictionary like object with parameter key/value pairs.

    Note that this class does not yet fully support reading all features
    of CABO files. For example, the parsing of booleans, date/times and
    tabular parameters is not supported and will lead to errors.

    The header of the CABO file (marked with ** at the first line) is
    read and can be retrieved by the get_header() method or just by
    a print on the returned dictionary.
    
    *Example*
    
    A parameter file 'parfile.cab' which looks like this::
    
        ** CROP DATA FILE for use with WOFOST Version 5.4, June 1992
        **
        ** WHEAT, WINTER 102
        ** Regions: Ireland, central en southern UK (R72-R79), 
        **          Netherlands (not R47), northern Germany (R11-R14)
        CRPNAM='Winter wheat 102, Ireland, N-U.K., Netherlands, N-Germany'
        CROP_NO=99
        TBASEM   = -10.0    ! lower threshold temp. for emergence [cel]
        DTSMTB   =   0.00,    0.00,     ! daily increase in temp. sum 
                    30.00,   30.00,     ! as function of av. temp. [cel; cel d]
                    45.00,   30.00
        ** maximum and minimum concentrations of N, P, and K
        ** in storage organs        in vegetative organs [kg kg-1]
        NMINSO   =   0.0110 ;       NMINVE   =   0.0030
    
    Can be read with the following statements::
    
        >>>fileparameters = CABOFileReader('parfile.cab')
        >>>print fileparameters['CROP_NO']
        99
        >>>print fileparameters
        ** CROP DATA FILE for use with WOFOST Version 5.4, June 1992
        **
        ** WHEAT, WINTER 102
        ** Regions: Ireland, central en southern UK (R72-R79),
        **          Netherlands (not R47), northern Germany (R11-R14)
        ------------------------------------
        TBASEM: -10.0 <type 'float'>
        DTSMTB: [0.0, 0.0, 30.0, 30.0, 45.0, 30.0] <type 'list'>
        NMINVE: 0.003 <type 'float'>
        CROP_NO: 99 <type 'int'>
        CRPNAM: Winter wheat 102, Ireland, N-U.K., Netherlands, N-Germany <type 'str'>
        NMINSO: 0.011 <type 'float'>
    """

    # RE patterns for parsing scalar, table and string parameters
    scpar = "[a-zA-Z0-9_]+[\s]*=[\s]*[a-zA-Z0-9_.\-]+"
    tbpar = "[a-zA-Z0-9_]+[\s]*=[\s]*[0-9,.\s\-+]+"
    strpar = "[a-zA-Z0-9_]+[\s]*=[\s]*'.*?'"

    def _remove_empty_lines(self, filecontents):
        t = []
        for line in filecontents:
            line = line.strip(" \n\r")
            if len(line)>0:
                t.append(line)
        return t
        
    def _remove_inline_comments(self, filecontents):
        t = []
        for line in filecontents:
            line = line.split("!")[0]
            line.strip()
            if len(line) > 0:
                t.append(line)
        return t
    
    def _is_comment(self, line):
        if line.startswith("*"):
            return True
        else:
            return False

    def _find_header(self, filecontents):
        """Parses and strips header marked with '*' at the beginning of 
        the file. Further lines marked with '*' are deleted."""

        header = []
        other_contents = []
        inheader = True
        for line in filecontents:
            if inheader is True:
                if self._is_comment(line):
                    header.append(line)
                else:
                    other_contents.append(line)
                    inheader = False
            else:
                if self._is_comment(line):
                    pass
                else:
                    other_contents.append(line)
        return (header, other_contents)
 
    def _parse_table_values(self, parstr):
        """Parses table parameter into a list of floats."""

        tmpstr = parstr.strip()
        valuestrs = tmpstr.split(",")
        if len(valuestrs) < 4:
            raise LengthError((len(valuestrs), valuestrs))
        if (len(valuestrs) % 2) != 0:
            raise XYPairsError((len(valuestrs), valuestrs))

        tblvalues = []
        for vstr in valuestrs:
            value = float(vstr)
            tblvalues.append(value)
        return tblvalues
        
    def _find_parameter_sections(self, filecontents):
        "returns the sections defining float, string and table parameters."
        scalars = ""
        strings = ""
        tables = ""
        
        for line in filecontents:
            if line.find("'") != -1: # string parameter
                strings += (line + " ")
            elif line.find(",") != -1: # table parameter
                tables += (line + " ")
            else:
                scalars += (line + " ")
                
        return scalars, strings, tables
       
    def _find_individual_pardefs(self, regexp, parsections):
        """Splits the string into individual parameters definitions.
        """
        par_definitions = re.findall(regexp, parsections)
        rest = re.sub(regexp, "", parsections)
        rest = rest.replace(";", "")
        if rest.strip() != "":
            msg = "Failed to parse the CABO file!\n" +\
                  ("Found the following parameter definitions:\n %s" % par_definitions) + \
                  ("Failed to parse:\n %s" % rest)
            raise PCSEError(msg)
        return par_definitions
        
    def __init__(self, fname):
        with open(fname) as fp:
            filecontents = fp.readlines()
        filecontents = self._remove_empty_lines(filecontents)
        filecontents = self._remove_inline_comments(filecontents)

        if len(filecontents) == 0:
            msg = "Empty CABO file!"
            raise PCSEError(msg)

        # Split between file header and parameters
        self.header, filecontents = self._find_header(filecontents)

        # Find parameter sections using string methods
        scalars, strings, tables = self._find_parameter_sections(filecontents)

        # Parse into individual parameter definitions
        scalar_defs = self._find_individual_pardefs(self.scpar, scalars)
        table_defs = self._find_individual_pardefs(self.tbpar, tables)
        string_defs = self._find_individual_pardefs(self.strpar, strings)

        # Parse individual parameter definitions into name & value.
        for parstr in scalar_defs:
            try:
                parname, valuestr = parstr.split("=")
                parname = parname.strip()
                if valuestr.find(".") != -1:
                    value = float(valuestr)
                else:
                    value = int(valuestr)
                self[parname] = value
            except (ValueError) as exc:
                msg = "Failed to parse parameter, value: %s, %s" 
                raise PCSEError(msg % (parstr, valuestr))

        for parstr in string_defs:
            try:
                parname, valuestr = parstr.split("=", 1)
                parname = parname.strip()
                value = (valuestr.replace("'","")).replace('"','')
                self[parname] = value
            except (ValueError) as exc:
                msg = "Failed to parse parameter, value: %s, %s" 
                raise PCSEError(msg % (parstr, valuestr))

        for parstr in table_defs:
            parname, valuestr = parstr.split("=")
            parname = parname.strip()
            try:
                value = self._parse_table_values(valuestr)
                self[parname] = value
            except (ValueError) as exc:
                msg = "Failed to parse table parameter %s: %s" % (parname, valuestr)
                raise PCSEError(msg)
            except (LengthError) as exc:
                msg = "Failed to parse table parameter %s: %s. \n" % (parname, valuestr)
                msg += "Table parameter should contain at least 4 values "
                msg += "instead got %i" 
                raise PCSEError(msg % exc.value[0])
            except (XYPairsError) as exc:
                msg = "Failed to parse table parameter %s: %s\n" % (parname, valuestr)
                msg += "Parameter should be have even number of positions."
                raise XYPairsError(msg)

    def __str__(self):
        msg = ""
        for line in self.header:
            msg += line+"\n"
        msg += "------------------------------------\n"
        for key, value in self.items():
            msg += ("%s: %s %s\n" % (key, value, type(value)))
        return msg