# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import SimulationObject, ParamTemplate, StatesTemplate, RatesTemplate
from pcse.crop.lintul_cassava.phenology import phenology
from pcse.crop.lintul_cassava.canopy_rain_interception import canopy_rain_interception
from pcse.crop.lintul_cassava.biomass_partitioning import biomass_partitioning
from pcse.crop.lintul_cassava.dormancy_and_recovery import dormancy_and_recovery
from pcse.crop.lintul_cassava.lintul_cassava_penman import penman
from pcse.crop.lintul_cassava.evapotranspiration import evapotranspiration
from pcse.crop.lintul_cassava.fibrous_root_growth import fibrous_root_growth
from pcse.crop.lintul_cassava.green_leaf_area import green_leaf_area
from pcse.crop.lintul_cassava.leaf_senescence import leaf_senescence
from pcse.crop.lintul_cassava.light_interception_and_growth import light_interception_and_growth
from pcse.crop.lintul_cassava.nutrient_dynamics import crop_nutrient_dynamics
from pcse.crop.lintul_cassava.nutrient_stress import npk_stress
from pcse.traitlets import Instance, Float

class LINTUL_CASSAVA(SimulationObject):
    """
    Top level object organizing the different components of the LINTUL Cassava
    simulation.

    The CropSimulation object organizes the different processes of the crop
    simulation. Moreover, it contains the parameters, rate and state variables
    which are relevant at the level of the entire crop. The processes that are
    implemented as embedded simulation objects consist of:

    1. Phenology (self.phenology)
    2. Fibrous root growth (self.fibrous_root_growth)
    3. Canopy rain interception (self.canopy_rain_interception)
    4. Evapotranspiration (self.evapotranspiration)
    5. NPK stress effects (self.npk_stress)
    6. Dormancy and recovery (self.dormancy)
    7. Leaf senescence (self.leaf_senescence)
    8. Light interception and growth (self.light_interception_and_growth)
    9. Biomass partitioning (self.biomass_partitioning)
    10. Crop nutrient dynamics (self.crop_nutrient_dynamics)
    11. Green leaf area dynamics (self.green_leaf_area)

    The LINTUL Cassava crop growth model in PCSE is a Python implementation of the LINTUL Cassava
    model (Ezui et al., 2018; Adiele et al., 2022) that was originally implemented in R.

    References:
    Adiele J.G., Schut A.G.T., Ezui K.S., Giller K.E. (2022) LINTUL-Cassava-NPK: A simulation
    model for nutrient-limited cassava growth. Field Crops Research 281: ARTN 108488.
    https://doi.org/10.1007/s13593-020-00649-w

    Ezui K.S., Leffelaar P.A., Franke A.C., Mando A., Giller K.E. (2018) Simulating drought impact
    and mitigation in cassava using the LINTUL model. Field Crops Research 219: 256-272.
    https://doi.org/10.1016/j.fcr.2018.01.033
    """
    phenology = Instance(SimulationObject)
    fibrous_root_growth = Instance(SimulationObject)
    penman = Instance(SimulationObject)
    canopy_rain_interception = Instance(SimulationObject)
    evapotranspiration = Instance(SimulationObject)
    npk_stress = Instance(SimulationObject)
    dormancy = Instance(SimulationObject)
    leaf_senescence = Instance(SimulationObject)
    light_interception_and_growth = Instance(SimulationObject)
    biomass_partitioning = Instance(SimulationObject)
    crop_nutrient_dynamics = Instance(SimulationObject)
    green_leaf_area = Instance(SimulationObject)

    class Parameters(ParamTemplate):
        pass

    class StateVariables(StatesTemplate):
        TAGP = Float()
        CTRAT = Float()
        CEVST = Float()

    class RateVariables(RatesTemplate):
        pass

    def initialize(self, day, kiosk, parvalues):
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.states = self.StateVariables(kiosk, TAGP=0.0, CTRAT=0.0, CEVST=0.0)
        self.rates = self.RateVariables(kiosk)

        self.phenology = phenology(day, kiosk, parvalues)
        self.fibrous_root_growth = fibrous_root_growth(day, kiosk, parvalues)
        self.canopy_rain_interception = canopy_rain_interception(day, kiosk, parvalues)
        self.evapotranspiration = evapotranspiration(day, kiosk, parvalues)
        self.npk_stress = npk_stress(day, kiosk, parvalues)
        self.dormancy = dormancy_and_recovery(day, kiosk, parvalues)
        self.leaf_senescence = leaf_senescence(day, kiosk, parvalues)
        self.light_interception_and_growth = light_interception_and_growth(day, kiosk, parvalues)
        self.biomass_partitioning = biomass_partitioning(day, kiosk, parvalues)
        self.crop_nutrient_dynamics = crop_nutrient_dynamics (day, kiosk, parvalues)
        self.green_leaf_area = green_leaf_area(day, kiosk, parvalues)

    def calc_rates(self, day, drv, delt = 1):
        self.phenology.calc_rates(day, drv, delt)
        self.fibrous_root_growth.calc_rates(day, drv, delt)
        self.canopy_rain_interception.calc_rates(day, drv, delt)
        self.evapotranspiration(day, drv)
        self.npk_stress(day, drv)
        self.dormancy.calc_rates(day, drv, delt)
        self.light_interception_and_growth.calc_rates(day, drv, delt)
        self.leaf_senescence.calc_rates(day, drv)
        self.biomass_partitioning.calc_rates(day, drv)
        self.crop_nutrient_dynamics.calc_rates(day, drv)
        self.green_leaf_area.calc_rates(day, drv)

    def integrate(self, day, delt = 1):
        self.phenology.integrate(day, delt)
        self.fibrous_root_growth.integrate(day, delt)
        self.canopy_rain_interception.integrate(day, delt)
        self.dormancy.integrate(day, delt)
        self.leaf_senescence.integrate(day, delt)
        self.light_interception_and_growth.integrate(day, delt)
        self.biomass_partitioning.integrate(day, delt)
        self.crop_nutrient_dynamics.integrate(day, delt)
        self.green_leaf_area.integrate(day, delt)

        k = self.kiosk
        s = self.states

        # Total above-ground biomass
        s.TAGP = k.WST + k.WLV + k.WSO

        # total crop transpiration and soil evaporation
        s.CTRAT += k.TRA * delt
        s.CEVST += k.EVS * delt


class LINTUL_CASSAVA_NO_NUTRIENT_STRESS(LINTUL_CASSAVA):
    """
    Top level object organizing the different components of the LINTUL Cassava
    simulation, assuming non-nutrient limited conditions.

    This class inherits all functionality from the LINTUL_CASSAVA class, except
    that the class variable NUTRIENT_LIMITED is set to False to prevent nutrient
    stress.

    The CropSimulation object organizes the different processes of the crop
    simulation. Moreover, it contains the parameters, rate and state variables
    which are relevant at the level of the entire crop. The processes that are
    implemented as embedded simulation objects consist of:

    1. Phenology (self.phenology)
    2. Fibrous root growth (self.fibrous_root_growth)
    3. Canopy rain interception (self.canopy_rain_interception)
    4. Evapotranspiration (self.evapotranspiration)
    5. NPK stress effects (self.npk_stress)
    6. Dormancy and recovery (self.dormancy)
    7. Leaf senescence (self.leaf_senescence)
    8. Light interception and growth (self.light_interception_and_growth)
    9. Biomass partitioning (self.biomass_partitioning)
    10. Crop nutrient dynamics (self.crop_nutrient_dynamics)
    11. Green leaf area dynamics (self.green_leaf_area)

    The LINTUL Cassava crop growth model in PCSE is a Python implementation of the LINTUL Cassava
    model (Ezui et al., 2018; Adiele et al., 2022) that was originally implemented in R.

    References:
    Adiele J.G., Schut A.G.T., Ezui K.S., Giller K.E. (2022) LINTUL-Cassava-NPK: A simulation
    model for nutrient-limited cassava growth. Field Crops Research 281: ARTN 108488.
    https://doi.org/10.1007/s13593-020-00649-w

    Ezui K.S., Leffelaar P.A., Franke A.C., Mando A., Giller K.E. (2018) Simulating drought impact
    and mitigation in cassava using the LINTUL model. Field Crops Research 219: 256-272.
    https://doi.org/10.1016/j.fcr.2018.01.033
    """
    def initialize(self, day, kiosk, parvalues):
        super().initialize(day, kiosk, parvalues)
        self.npk_stress.NUTRIENT_LIMITED = False

class LINTUL_CASSAVA_original(LINTUL_CASSAVA):
    """
    Top level object organizing the different components of the LINTUL Cassava
    simulation, assuming non-nutrient limited conditions.

    This class inherits all functionality from the LINTUL_CASSAVA class, except
    for that it calculates the reference evapotranspiration rate with its native
    Penman method, instead of using the PCSE functionalities to calculate this.

    The CropSimulation object organizes the different processes of the crop
    simulation. Moreover, it contains the parameters, rate and state variables
    which are relevant at the level of the entire crop. The processes that are
    implemented as embedded simulation objects consist of:

    1. Phenology (self.phenology)
    2. Fibrous root growth (self.fibrous_root_growth)
    3. Canopy rain interception (self.canopy_rain_interception)
    4. Evapotranspiration (self.evapotranspiration)
    5. NPK stress effects (self.npk_stress)
    6. Dormancy and recovery (self.dormancy)
    7. Leaf senescence (self.leaf_senescence)
    8. Light interception and growth (self.light_interception_and_growth)
    9. Biomass partitioning (self.biomass_partitioning)
    10. Crop nutrient dynamics (self.crop_nutrient_dynamics)
    11. Green leaf area dynamics (self.green_leaf_area)

    The LINTUL Cassava crop growth model in PCSE is a Python implementation of the LINTUL Cassava
    model (Ezui et al., 2018; Adiele et al., 2022) that was originally implemented in R.

    References:
    Adiele J.G., Schut A.G.T., Ezui K.S., Giller K.E. (2022) LINTUL-Cassava-NPK: A simulation
    model for nutrient-limited cassava growth. Field Crops Research 281: ARTN 108488.
    https://doi.org/10.1007/s13593-020-00649-w

    Ezui K.S., Leffelaar P.A., Franke A.C., Mando A., Giller K.E. (2018) Simulating drought impact
    and mitigation in cassava using the LINTUL model. Field Crops Research 219: 256-272.
    https://doi.org/10.1016/j.fcr.2018.01.033
    """
    def initialize(self, day, kiosk, parvalues):
        super().initialize(day, kiosk, parvalues)
        self.penman = penman(day, kiosk, parvalues)

    def calc_rates(self, day, drv, delt = 1):
        self.penman(day, drv)
        k = self.kiosk
        drv.ET0 = k.ET0
        drv.ES0 = k.ES0
        super().calc_rates(day, drv, delt = 1)

class LINTUL_CASSAVA_original_NO_NUTRIENT_STRESS(LINTUL_CASSAVA_original):
    """
    Top level object organizing the different components of the LINTUL Cassava
    simulation, assuming non-nutrient limited conditions.

    This class inherits all functionality from the LINTUL_CASSAVA class, except
    for that it calculates the reference evapotranspiration rate with its native
    Penman method, instead of using the PCSE functionalities to calcuate this and
    that it assumes no-nutrient limited conditions. The latter is obtained by setting
    the class variable NUTRIENT_LIMITED to false.

    The CropSimulation object organizes the different processes of the crop
    simulation. Moreover, it contains the parameters, rate and state variables
    which are relevant at the level of the entire crop. The processes that are
    implemented as embedded simulation objects consist of:

    1. Phenology (self.phenology)
    2. Fibrous root growth (self.fibrous_root_growth)
    3. Canopy rain interception (self.canopy_rain_interception)
    4. Evapotranspiration (self.evapotranspiration)
    5. NPK stress effects (self.npk_stress)
    6. Dormancy and recovery (self.dormancy)
    7. Leaf senescence (self.leaf_senescence)
    8. Light interception and growth (self.light_interception_and_growth)
    9. Biomass partitioning (self.biomass_partitioning)
    10. Crop nutrient dynamics (self.crop_nutrient_dynamics)
    11. Green leaf area dynamics (self.green_leaf_area)

    The LINTUL Cassava crop growth model in PCSE is a Python implementation of the LINTUL Cassava
    model (Ezui et al., 2018; Adiele et al., 2022) that was originally implemented in R.

    References:
    Adiele J.G., Schut A.G.T., Ezui K.S., Giller K.E. (2022) LINTUL-Cassava-NPK: A simulation
    model for nutrient-limited cassava growth. Field Crops Research 281: ARTN 108488.
    https://doi.org/10.1007/s13593-020-00649-w

    Ezui K.S., Leffelaar P.A., Franke A.C., Mando A., Giller K.E. (2018) Simulating drought impact
    and mitigation in cassava using the LINTUL model. Field Crops Research 219: 256-272.
    https://doi.org/10.1016/j.fcr.2018.01.033
    """
    def initialize(self, day, kiosk, parvalues):
        super().initialize(day, kiosk, parvalues)
        self.npk_stress.NUTRIENT_LIMITED = False