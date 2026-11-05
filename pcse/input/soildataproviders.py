import copy
import yaml


class DummySoilDataProvider(dict):
    """This class is to provide some dummy soil parameters for potential production simulation.

    Simulation of potential production levels is independent of the soil. Nevertheless, the model
    does not some parameter values. This data provider provides some hard coded parameter values for
    this situation.
    """
    _defaults = {"SMFCF":0.3,
                 "SM0":0.4,
                 "SMW":0.1,
                 "RDMSOL":120,
                 "CRAIRC":0.06,
                 "K0":10.,
                 "SOPE":10.,
                 "KSUB":10.}

    def __init__(self):
        dict.__init__(self)
        self.update(self._defaults)

    def copy(self):
        """
        Overrides the inherited dict.copy method, which returns a dict.
        This instead preserves the class and attributes like .header.
        """
        return copy.copy(self)


class CWBSimpleSoilDataProvider(dict):
    """This provides soil parameter definitions for some generic soil types for the classic water balance
    derived from the EU soil map.

    This is simply to have some default parameters for soil types for demonstration purposes.

    Starting the CWBSimpleSoilDataProvider with no arguments prints the available soil names.

    :param soil_name: the name of the soil
    :param rootable_depth: the rootable depth of the soil in cm
    """
    soil_definitions = yaml.safe_load("""
        default_parameters: &defaults
            K0:   10.0   # hydraulic conductivity of saturated soil [cm day-1]
            SOPE: 10.0   # maximum percolation rate root zone[cm day-1]
            KSUB: 10.0   # maximum percolation rate subsoil [cm day-1]
        texture_1_coarse:
            <<: *defaults
            SMW:    0.040  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF:  0.190  #  soil moisture content at field capacity [cm3/cm3]
            SM0:    0.390  #  soil moisture content at saturation [cm3/cm3]           
            CRAIRC: 0.090  # critical soil air content for aeration [cm3/cm3]
        texture_2_medium_coarse:
            <<: *defaults
            SMW:   0.100  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF: 0.320  #  soil moisture content at field capacity [cm3/cm3]
            SM0:   0.430  #  soil moisture content at saturation [cm3/cm3]
            CRAIRC: 0.075 # critical soil air content for aeration [cm3 /cm3]
        texture_3_medium_fine:
            <<: *defaults
            SMW:   0.100  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF: 0.350  #  soil moisture content at field capacity [cm3/cm3]
            SM0:   0.450  #  soil moisture content at saturation [cm3/cm3]
            CRAIRC: 0.075 # critical soil air content for aeration [cm3 /cm3]
        texture_4_fine:
            <<: *defaults
            SMW:   0.160  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF: 0.350  #  soil moisture content at field capacity [cm3/cm3]
            SM0:   0.450  #  soil moisture content at saturation [cm3/cm3]
            CRAIRC: 0.075 # critical soil air content for aeration [cm3 /cm3]
        texture_5_very fine:
            <<: *defaults
            SMW:   0.220  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF: 0.350  #  soil moisture content at field capacity [cm3/cm3]
            SM0:   0.450  #  soil moisture content at saturation [cm3/cm3]
            CRAIRC: 0.075 # critical soil air content for aeration [cm3 /cm3]
        texture_6_peat:
            <<: *defaults
            SMW:    0.130  #  soil moisture content at wilting point [cm3/cm3]
            SMFCF:  0.522  #  soil moisture content at field capacity [cm3/cm3]
            SM0:    0.730  #  soil moisture content at saturation [cm3/cm3]
            CRAIRC: 0.100  # critical soil air content for aeration [cm3/cm3]
        """)

    def __init__(self, soil_name=None, rootable_depth=120.0):
        super().__init__()
        soils = [v for v in self.soil_definitions.keys() if not v.startswith("default")]
        if soil_name is None:
            print(f"Error: Soil name not provide, should be provided as one of: {soils}")
            return
        elif soil_name not in soils:
            print(f"Error: Unrecognized soil: '{soil_name}'. Should be one of: {soils}")
            return

        self.update(self.soil_definitions[soil_name])

        if 20. <= rootable_depth <= 250.:
            self.update({"RDMSOL": rootable_depth})
        else:
            print(f"Error: Rootable depth ({rootable_depth}) must be between 20 and 250 cm.")



