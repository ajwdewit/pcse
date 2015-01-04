# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from ..traitlets import Float, Int, Instance, Dict, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, SimulationObject

class WOFOST_Maintenance_Respiration(SimulationObject):
    """Maintenance respiration in WOFOST
    
    WOFOST calculates the maintenance respiration as proportional to the dry
    weights of the plant organs to be maintained, where each plant organ can be
    assigned a different maintenance coefficient. Multiplying organ weight
    with the maintenance coeffients yields the relative maintenance respiration
    (`RMRES`) which is than corrected for senescence (parameter `RFSETB`). Finally,
    the actual maintenance respiration rate is calculated using the daily mean
    temperature, assuming a relative increase for each 10 degrees increase
    in temperature as defined by `Q10`.

    **Simulation parameters:** (To be provided in cropdata dictionary):
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    Q10      Relative increase in maintenance repiration    SCr       -
             rate with each 10 degrees increase in
             temperature
    RMR      Relative maintenance respiration rate for
             roots                                          SCr     |kg CH2O kg-1 d-1|
    RMS      Relative maintenance respiration rate for
             stems                                          SCr     |kg CH2O kg-1 d-1|
    RML      Relative maintenance respiration rate for
             leaves                                         SCr     |kg CH2O kg-1 d-1|
    RMO      Relative maintenance respiration rate for
             storage organs                                 SCr     |kg CH2O kg-1 d-1|
    =======  ============================================= =======  ============
    

    **State and rate variables:**
    
    `WOFOSTMaintenanceRespiration` has no state/rate variables, but calculates
    the rate of respiration which is returned directly from the `__call__()`
    method.
    
    **Signals send or handled**
    
    None
    
    **External dependencies:**
    
    =======  =================================== =============================  ============
     Name     Description                         Provided by                    Unit
    =======  =================================== =============================  ============
    DVS      Crop development stage              DVS_Phenology                  -
    WRT      Dry weight of living roots          WOFOST_Root_Dynamics           |kg ha-1|
    WST      Dry weight of living stems          WOFOST_Stem_Dynamics           |kg ha-1|
    WLV      Dry weight of living leaves         WOFOST_Leaf_Dynamics           |kg ha-1|
    WSO      Dry weight of living storage organs WOFOST_Storage_Organ_Dynamics  |kg ha-1|
    =======  =================================== =============================  ============


    """
    
    class Parameters(ParamTemplate):
        Q10 = Float(-99.)
        RMR = Float(-99.)
        RML = Float(-99.)
        RMS = Float(-99.)
        RMO = Float(-99.)
        RFSETB = AfgenTrait()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk
        
    def __call__(self, day, drv):
        p = self.params
        kk = self.kiosk
        
        RMRES = (p.RMR * kk["WRT"] + \
                 p.RML * kk["WLV"] + \
                 p.RMS * kk["WST"] + \
                 p.RMO * kk["WSO"])
        RMRES *= p.RFSETB(kk["DVS"])
        TEFF   = p.Q10**((drv.TEMP-25.)/10.)
        return RMRES*TEFF
