# Define single and tabular crop parameter values
WOFOST_parameter_codes_single = ("CFET", "CVL", "CVO", "CVR", "CVS", "DEPNR", "DLC",
                          "DLO", "DVSEND", "EFF", "IAIRDU", "IDSL", "KDIF",
                          "LAIEM", "PERDL", "Q10", "RDI", "RDMCR", "RGRLAI",
                          "RML", "RMO", "RMR", "RMS", "RRI", "SPA", "SPAN", "SSA",
                          "TBASE", "TBASEM", "TDWI", "TEFFMX", "TSUM1", "TSUM2",
                          "TSUMEM", "VERNSAT", "VERNBASE", "VERNDVS")
WOFOST_parameter_codes_tabular = ("AMAXTB", "DTSMTB", "FLTB", "FOTB", "FRTB", "FSTB",
                           "RDRRTB", "RDRSTB", "RFSETB", "SLATB", "TMNFTB",
                           "TMPFTB", "VERNRTB")
# Some parameters have to be converted from a single to a tabular form
WOFOST_single2tabular = {"SSA": ("SSATB", [0., None, 2.0, None]),
                         "KDIF": ("KDIFTB", [0., None, 2.0, None]),
                         "EFF": ("EFFTB", [0., None, 40., None])}
# Default values for additional parameters not defined in CGMS
WOFOST_parameters_additional = {"DVSI": 0.0, "IOX": 0}
# Parameters which are optional (mainly dealing with vernalisation)
WOFOST_optional_parameters = ["VERNSAT", "VERNBASE", "VERNDVS", "VERNRTB"]
