# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
from .. import exceptions as exc


class _GenericSiteDataProvider(dict):
    """Generic Site data provider

    It just checks the values provided as keywords given the _defaults and _required values

    _defaults = {"VARNAME": (default, {maxvalue, minvalue}, type),
                }
    _required = ["VARNAME"]
    """

    def __init__(self, **kwargs):
        dict.__init__(self)

        for par_name, (default_value, par_range, par_conversion) in self._defaults.items():
            if par_name not in kwargs:
                # parameter was not provided, use the default if possible
                if par_name in self._required:
                    msg = "Value for parameter '%s' must be provided!" % par_name
                    raise exc.PCSEError(msg)
                else:
                    par_value = default_value
            else:
                # parameter was provided, check value for type and range
                par_value = par_conversion(kwargs.pop(par_name))
                if isinstance(par_range, set):
                    # allowed values consist of a
                    if par_value not in par_range:
                        msg = "Value for parameter '%s' can only have values: %s" % (par_name, par_range)
                        raise exc.PCSEError(msg)
                else:
                    if isinstance(par_value, list):
                        if not all(par_range[0] <= x <= par_range[1] for x in par_value):
                            msg = "At least one of the values for parameter '%s' out of range %s-%s" % \
                                  (par_name, par_range[0], par_range[1])
                            raise exc.PCSEError(msg)
                    else:
                        if not (par_range[0] <= par_value <= par_range[1]):
                            msg = "Value for parameter '%s' out of range %s-%s" % \
                                  (par_name, par_range[0], par_range[1])
                            raise exc.PCSEError(msg)
            self[par_name] = par_value

        # Check if kwargs is empty
        if kwargs:
            msg = f"Unknown parameter values provided to {self.__class__}: %s" % kwargs
            raise exc.PCSEError(msg)


class WOFOST72SiteDataProvider(_GenericSiteDataProvider):
    """Site data provider for WOFOST 7.2.

    Site specific parameters for WOFOST 7.2 can be provided through this data provider as well as through
    a normal python dictionary. The sole purpose of implementing this data provider is that the site
    parameters for WOFOST are documented, checked and that sensible default values are given.

    The following site specific parameter values can be set through this data provider::

        - IFUNRN    Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
                    or not (0). Default 0
        - NOTINF    Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
        - SSMAX     Maximum depth of water that can be stored on the soil surface [cm]
        - SSI       Initial depth of water stored on the surface [cm]
        - WAV       Initial amount of water in total soil profile [cm]
        - SMLIM     Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4

    Currently only the value for WAV is mandatory to specify.
    """

    _defaults = {"IFUNRN": (0, {0, 1}, int),
                 "NOTINF": (0, [0., 1.], float),
                 "SSI": (0., [0., 100.], float),
                 "SSMAX": (0., [0., 100.], float),
                 "WAV": (None, [0., 100.], float),
                 "SMLIM": (0.4, [0., 1.], float)}
    _required = ["WAV"]


class WOFOST73SiteDataProvider(_GenericSiteDataProvider):
    """Site data provider for WOFOST 7.3

    Site specific parameters for WOFOST 7.3 can be provided through this data provider as well as through
    a normal python dictionary. The sole purpose of implementing this data provider is that the site
    parameters for WOFOST are documented, checked and that sensible default values are given.

    The following site specific parameter values can be set through this data provider::

        - IFUNRN    Indicates whether non-infiltrating fraction of rain is a function of storm size (1)
                    or not (0). Default 0
        - NOTINF    Maximum fraction of rain not-infiltrating into the soil [0-1], default 0.
        - SSMAX     Maximum depth of water that can be stored on the soil surface [cm]
        - SSI       Initial depth of water stored on the surface [cm]
        - WAV       Initial amount of water in total soil profile [cm]
        - SMLIM     Initial maximum moisture content in initial rooting depth zone [0-1], default 0.4
        - CO2       Atmospheric CO2 concentration in ppm

    Values for WAV and CO2 is mandatory to specify.
    """

    _defaults = {"IFUNRN": (0, {0, 1}, int),
                 "NOTINF": (0, [0., 1.], float),
                 "SSI": (0., [0., 100.], float),
                 "SSMAX": (0., [0., 100.], float),
                 "WAV": (None, [0., 100.], float),
                 "SMLIM": (0.4, [0., 1.], float),
                 "CO2": (None, [320, 700], float)}
    _required = ["WAV", "CO2"]


class WOFOST81SiteDataProvider_Classic(_GenericSiteDataProvider):
    """Site data provider for WOFOST 8.1 for Classic water and nitrogen balance.

    Site specific parameters for WOFOST 8.1 can be provided through this data provider as well as through
    a normal python dictionary. The sole purpose of implementing this data provider is that the site
    parameters for WOFOST are documented, checked and that sensible default values are given.

    The following site specific parameter values can be set through this data provider::

        - IFUNRN        Indicates whether non-infiltrating fraction of rain is a function of
                        storm size (1) or not (0). Default 0
        - NOTINF        Maximum fraction of rain not-infiltrating into the soil [0-1],
                        default 0.
        - SSMAX         Maximum depth of water that can be stored on the soil surface [cm]
        - SSI           Initial depth of water stored on the surface [cm]
        - WAV           Initial amount of water in total soil profile [cm]
        - SMLIM         Initial maximum moisture content in initial rooting depth zone [0-1],
                        default 0.4
        - CO2           Atmospheric CO2 level (ppm), default 360.
        - BG_N_SUPPLY   Background N supply through atmospheric deposition in kg/ha/day. Can be
                        in the order of 25 kg/ha/year in areas with high N pollution. Default 0.0
        - NSOILBASE     Base N amount available in the soil. This is often estimated as the nutrient
                        left over from the previous growth cycle (surplus nutrients, crop residues
                        or green manure).
        - NSOILBASE_FR  Daily fraction of soil N coming available through mineralization
        - NAVAILI       Amount of N available in the pool at initialization of the system [kg/ha]

    Currently, the parameters for initial water availability (WAV) and initial availability of
    nutrients (NAVAILI) are mandatory to specify.
    """

    _defaults = {"IFUNRN": (0, {0, 1}, int),
                 "NOTINF": (0, [0., 1.], float),
                 "SSI": (0., [0., 100.], float),
                 "SSMAX": (0., [0., 100.], float),
                 "WAV": (None, [0., 100.], float),
                 "SMLIM": (0.4, [0., 1.], float),
                 "CO2": (None, [300., 1400.], float),
                 "BG_N_SUPPLY": (0, (0, 0.1), float),
                 "NSOILBASE": (0, (0, 100), float),
                 "NSOILBASE_FR": (0.025, (0, 100), float),
                 "NAVAILI": (None, (0, 250), float),
                 }
    _required = ["WAV", "NAVAILI", "CO2"]


class WOFOST81SiteDataProvider_SNOMIN(_GenericSiteDataProvider):
    """Site data provider for WOFOST 8.1 for use with the SNOMIN C/N balance.

    The following site specific parameter values can be set through this data provider::

        - IFUNRN        Indicates whether non-infiltrating fraction of rain is a function of
                        storm size (1) or not (0). Default 0
        - NOTINF        Maximum fraction of rain not-infiltrating into the soil [0-1],
                        default 0.
        - SSMAX         Maximum depth of water that can be stored on the soil surface [cm]
        - SSI           Initial depth of water stored on the surface [cm]
        - WAV           Initial amount of water in total soil profile [cm]
        - SMLIM         Initial maximum moisture content in initial rooting depth zone [0-1],
                        default 0.4
        - CO2           Atmospheric CO2 level, currently around 400. [ppm]
        - A0SOM         Initial age of organic material (24.0)  [year]
        - CNRatioBio    C:N ratio of microbial biomass  (9.0) [kg C kg-1 N]
        - FASDIS        Assimilation to dissimilation rate ratio (0.5) [-]
        - KDENIT_REF    Reference first order rate of denitrification (0.06) [d-1]
        - KNIT_REF      Reference first order rate of nitrification (1.0) [d-1]
        - KSORP         Sorption coefficient (0.0005) [m3 soil kg-1 soil]
        - MRCDIS        Michaelis-Menten constant of relationship organic C-dissimilation rate
                        and response factor denitrification rate (0.001) [kg C m-2 d-1]
        - NH4ConcR      NH4-N concentration in rain water (0.9095) [mg NH4+-N L-1 water]
        - NO3ConcR      NO3-N concentration in rain water (2.1) [mg NO3--N L-1 water]
        - NH4I          Initial amount of NH4+ per soil layer  [kg NH4+ ha-1]. This
                        should match the number of soil layers specified in the soil
                        configuration. The initial value can be highly variable and as
                        high as 300-500 kg/ha of NH4/NO3 if the model was started right
                        after an N application event.
        - NO3I          Initial amount of NO3-N per soil layer [kg NO3-N ha-1]. This
                        should match the number of soil layers specified in the soil
                        configuration. The initial value can be highly variable and as
                        high as 300-500 kg/ha of NH4/NO3 if the model was started right
                        after an N application event.
        - WFPS_CRIT     Critical fraction water filled soil pores (0.8)  [m3 water m-3 pores]


        *important*: Some of the valid ranges of parameters for WOFOST 8.1/SNOMIN are uncertain
        and therefore values outside of the specified ranges here may be valid in certain cases.

    """
    _required = ["WAV", "CO2", "NH4I", "NO3I", ]
    _defaults = {"IFUNRN": (0, {0, 1}, int),
                 "NOTINF": (0, [0., 1.], float),
                 "SSI": (0., [0., 100.], float),
                 "SSMAX": (0., [0., 100.], float),
                 "WAV": (None, [0., 100.], float),
                 "SMLIM": (0.4, [0., 1.], float),
                 "CO2": (None, [300., 1400.], float),
                 "A0SOM": (24.0, [5.0, 40.0], float),
                 "CNRatioBio": (9.0, [5.0, 20.0], float),
                 "FASDIS": (0.5, [0, 0.6], float),
                 "KDENIT_REF": (0.06, [0.0, 0.1], float),
                 "KNIT_REF": (1.0, [0.9, 1.0], float),
                 "KSORP": (0.0005, [0.0001, 0.001], float),
                 "MRCDIS": (0.001, [0.0001, 0.01], float),
                 "NH4ConcR": (0.0, [0.0, 5.], float),
                 "NO3ConcR": (0.0, [0.0, 20.], float),
                 "NH4I": (None, [0.0, 300.0], list),
                 "NO3I": (None, [0.0, 500.0], list),
                 "WFPS_CRIT": (0.8, [0.5, 0.99], float),
                 }
