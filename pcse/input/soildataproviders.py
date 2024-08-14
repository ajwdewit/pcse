import copy

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
