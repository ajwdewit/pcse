# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
#
# Adding the impact of N on partitioning
# Herman Berghuijs (herman.berghuijs@wur.nl), April 2023
from collections import namedtuple
from math import exp

from ..traitlets import Float, Int, Instance
from ..decorators import prepare_rates, prepare_states
from ..base import ParamTemplate, StatesTemplate, SimulationObject,\
     VariableKiosk
from .. import exceptions as exc
from warnings import warn
from ..util import AfgenTrait


# Template for namedtuple containing partitioning factors
class PartioningFactors(namedtuple("partitioning_factors", "FR FL FS FO")):
    pass

class DVS_Partitioning(SimulationObject):
    """Class for assimilate partioning based on development stage (`DVS`).

    `DVS_partioning` calculates the partitioning of the assimilates to roots,
    stems, leaves and storage organs using fixed partitioning tables as a
    function of crop development stage. The available assimilates are first
    split into below-ground and abovegrond using the values in FRTB. In a
    second stage they are split into leaves (`FLTB`), stems (`FSTB`) and storage
    organs (`FOTB`).
    
    Since the partitioning fractions are derived from the state variable `DVS`
    they are regarded state variables as well.

    **Simulation parameters** (To be provided in cropdata dictionary):
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    FRTB     Partitioning to roots as a function of          TCr       -
             development stage.
    FSTB     Partitioning to stems as a function of          TCr       -
             development stage.
    FLTB     Partitioning to leaves as a function of         TCr       -
             development stage.
    FOTB     Partitioning to storage organs as a function    TCr       -
             of development stage.
    =======  ============================================= =======  ============
    

    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    FR        Fraction partitioned to roots.                     Y    -
    FS        Fraction partitioned to stems.                     Y    -
    FL        Fraction partitioned to leaves.                    Y    -
    FO        Fraction partitioned to storage orgains            Y    -
    =======  ================================================= ==== ============

    **Rate variables**

    None
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    =======  =================================== =================  ============
    
    *Exceptions raised*
    
    A PartitioningError is raised if the partitioning coefficients to leaves,
    stems and storage organs on a given day do not add up to '1'.
    """
    
    class Parameters(ParamTemplate):
        FRTB   = AfgenTrait()
        FLTB   = AfgenTrait()
        FSTB   = AfgenTrait()
        FOTB   = AfgenTrait()

    class StateVariables(StatesTemplate):
        FR = Float(-99.)
        FL = Float(-99.)
        FS = Float(-99.)
        FO = Float(-99.)
        PF = Instance(PartioningFactors)
    
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk

        # initial partitioning factors (pf)
        DVS = self.kiosk["DVS"]
        FR = self.params.FRTB(DVS)
        FL = self.params.FLTB(DVS)
        FS = self.params.FSTB(DVS)
        FO = self.params.FOTB(DVS)

        # Pack partitioning factors into tuple
        PF = PartioningFactors(FR, FL, FS, FO)
        
        # Initial states
        self.states = self.StateVariables(kiosk, publish=["FR","FL","FS","FO"],
                                          FR=FR, FL=FL, FS=FS, FO=FO, PF=PF)
        self._check_partitioning()
        
    def _check_partitioning(self):
        """Check for partitioning errors."""

        FR = self.states.FR
        FL = self.states.FL
        FS = self.states.FS
        FO = self.states.FO
        checksum = FR+(FL+FS+FO)*(1.-FR) - 1.
        if abs(checksum) >= 0.0001:
            msg = ("Error in partitioning!\n")
            msg += ("Checksum: %f, FR: %5.3f, FL: %5.3f, FS: %5.3f, FO: %5.3f\n" \
                    % (checksum, FR, FL, FS, FO))
            self.logger.error(msg)
            warn(msg)
#             raise exc.PartitioningError(msg)

    @prepare_states
    def integrate(self, day, delt=1.0):
        """Update partitioning factors based on development stage (DVS)"""

        params = self.params
        
        DVS = self.kiosk["DVS"]
        self.states.FR = params.FRTB(DVS)
        self.states.FL = params.FLTB(DVS)
        self.states.FS = params.FSTB(DVS)
        self.states.FO = params.FOTB(DVS)
        
        # Pack partitioning factors into tuple
        self.states.PF = PartioningFactors(self.states.FR, self.states.FL,
                                           self.states.FS, self.states.FO)

        self._check_partitioning()  
    
    def calc_rates(self, day, drv):
        """ Return partitioning factors based on current DVS.
        """
        # rate calculation does nothing for partioning as it is a derived
        # state
        return self.states.PF


# TODO: Ask Herman about this module, seems not doing anything with N stress
class DVS_Partitioning_N(SimulationObject):
    """Class for assimilate partitioning based on development stage (`DVS`)
    with influence of N stress.

    `DVS_Partitioning_N` calculates the partitioning of the assimilates to roots,
    stems, leaves and storage organs using fixed partitioning tables as a
    function of crop development stage. The only different with the normal
    partitioning class is the effect of nitrogen stress on partitioning to
    leaves (parameter NPART). The available assimilates are first
    split into below-ground and aboveground using the values in FRTB. In a
    second stage they are split into leaves (`FLTB`), stems (`FSTB`) and storage
    organs (`FOTB`).

    Since the partitioning fractions are derived from the state variable `DVS`
    they are regarded state variables as well.

    **Simulation parameters** (To be provided in cropdata dictionary):

    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    FRTB     Partitioning to roots as a function of          TCr       -
             development stage.
    FSTB     Partitioning to stems as a function of          TCr       -
             development stage.
    FLTB     Partitioning to leaves as a function of         TCr       -
             development stage.
    FOTB     Partitioning to starge organs as a function     TCr       -
             of development stage.
    NPART    Coefficient for the effect of N stress on       SCR       -
             leaf biomass allocation
    =======  ============================================= =======  ============


    **State variables**

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    FR        Fraction partitioned to roots.                     Y    -
    FS        Fraction partitioned to stems.                     Y    -
    FL        Fraction partitioned to leaves.                    Y    -
    FO        Fraction partitioned to storage orgains            Y    -
    =======  ================================================= ==== ============

    **Rate variables**

    None

    **Signals send or handled**

    None

    **External dependencies:**

    =======  =================================== =================  ============
     Name     Description                         Provided by         Unit
    =======  =================================== =================  ============
    DVS      Crop development stage              DVS_Phenology       -
    TRA      Actual transpiration                Simple_Evapotranspiration mm d-1
    TRAMX    Maximum transpiration               Simple_Evapotranspiration mm d-1
    NNI      Nitrogen nutrition index            npk_dynamics        -
    =======  =================================== =================  ============

    *Exceptions raised*

    A PartitioningError is raised if the partitioning coefficients to leaves,
    stems and storage organs on a given day do not add up to '1'.
    """

    class Parameters(ParamTemplate):
        FRTB = AfgenTrait()
        FLTB = AfgenTrait()
        FSTB = AfgenTrait()
        FOTB = AfgenTrait()
        #NPART = Float(-99.)  # coefficient for the effect of N stress on leaf allocation

    class StateVariables(StatesTemplate):
        FR = Float(-99.)
        FL = Float(-99.)
        FS = Float(-99.)
        FO = Float(-99.)
        PF = Instance(PartioningFactors)

    def initialize(self, day, kiosk, parameters):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parameters: dictionary with WOFOST cropdata key/value pairs
        """
        self.params = self.Parameters(parameters)

        # initial partioning factors (pf)
        k = self.kiosk
        FR = self.params.FRTB(k.DVS)
        FL = self.params.FLTB(k.DVS)
        FS = self.params.FSTB(k.DVS)
        FO = self.params.FOTB(k.DVS)

        # Pack partitioning factors into tuple
        PF = PartioningFactors(FR, FL, FS, FO)

        # Initial states
        self.states = self.StateVariables(kiosk, publish=["FR","FL","FS","FO"],
                                          FR=FR, FL=FL, FS=FS, FO=FO, PF=PF)
        self._check_partitioning()

    def _check_partitioning(self):
        """Check for partitioning errors."""

        FR = self.states.FR
        FL = self.states.FL
        FS = self.states.FS
        FO = self.states.FO
        checksum = FR+(FL+FS+FO)*(1.-FR) - 1.
        if abs(checksum) >= 0.0001:
            msg = ("Error in partitioning!\n")
            msg += ("Checksum: %f, FR: %5.3f, FL: %5.3f, FS: %5.3f, FO: %5.3f\n" \
                    % (checksum, FR, FL, FS, FO))
            self.logger.error(msg)
            raise exc.PartitioningError(msg)

    @prepare_states
    def integrate(self, day, delt=1.0):
        """
        Update partitioning factors based on development stage (DVS) and water and oxygen stress
        """

        p = self.params
        s = self.states
        k = self.kiosk

        FRTMOD = max(1., 1./(k.RFTRA + 0.5))
        s.FR = min(0.6, p.FRTB(k.DVS) * FRTMOD)
        s.FL = p.FLTB(k.DVS)
        s.FS = p.FSTB(k.DVS)
        s.FO = p.FOTB(k.DVS)

        # Pack partitioning factors into tuple
        s.PF = PartioningFactors(s.FR, s.FL, s.FS, s.FO)

        self._check_partitioning()

    def calc_rates(self, day, drv):
        """ Return partitioning factors based on current DVS.
        """
        # rate calculation does nothing for partitioning as it is a derived state
        return self.states.PF
