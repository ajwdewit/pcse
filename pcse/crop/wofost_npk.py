#!/usr/bin/env python

import datetime

from ..traitlets import Float, Int, Instance, Enum, Unicode
from ..decorators import prepare_rates, prepare_states
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals
from .. import exceptions as exc
from .phenology import DVS_Phenology as Phenology
from .respiration import WOFOST_Maintenance_Respiration as MaintenanceRespiration
from .stem_dynamics import WOFOST_Stem_Dynamics as Stem_Dynamics
from .root_dynamics import WOFOST_Root_Dynamics as Root_Dynamics
from .leaf_dynamics import WOFOST_Leaf_Dynamics_NPK as Leaf_Dynamics
from .storage_organ_dynamics import WOFOST_Storage_Organ_Dynamics as \
    Storage_Organ_Dynamics
from .assimilation import WOFOST_Assimilation2 as Assimilation
from .partitioning import DVS_Partitioning_NPK as Partitioning
from .evapotranspiration import Evapotranspiration2 as Evapotranspiration

from .npk_dynamics import NPK_Crop_Dynamics as NPK_crop
from .nutrients.npk_soil_dynamics import NPK_Soil_Dynamics as NPK_soil
from .nutrients.npk_stress import NPK_Stress as NPK_Stress


class WofostNPK(SimulationObject):
    
    """Top level object organizing the different components of the WOFOST crop
    simulation including the implementation of N/P/K dynamics.
            
    The CropSimulation object organizes the different processes of the crop
    simulation. Moreover, it contains the parameters, rate and state variables
    which are relevant at the level of the entire crop. The processes that are
    implemented as embedded simulation objects consist of:
    
        1. Phenology (self.pheno)
        2. Partitioning (self.part)
        3. Assimilation (self.assim)
        4. Maintenance respiration (self.mres)
        5. Evapotranspiration (self.evtra)
        6. Leaf dynamics (self.lv_dynamics)
        7. Stem dynamics (self.st_dynamics)
        8. Root dynamics (self.ro_dynamics)
        9. Storage organ dynamics (self.so_dynamics)
        10. N/P/K crop dynamics (self.npk_crop_dynamics)
        11. N/P/K soil dynamics (self.npk_soil_dynamics)
        12. N/P/K stress (self.npk_stress)

    **Simulation parameters:**
    
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    CVL      Conversion factor for assimilates to leaves       SCr     -
    CVO      Conversion factor for assimilates to storage      SCr     -
             organs.
    CVR      Conversion factor for assimilates to roots        SCr     -
    CVS      Conversion factor for assimilates to stems        SCr     -
    ======== =============================================== =======  ==========
    
    
    **State variables:**

    ============  ================================================= ==== ===============
     Name          Description                                      Pbl      Unit
    ============  ================================================= ==== ===============
    TAGP          Total above-ground Production                      N    |kg ha-1|
    GASST         Total gross assimilation                           N    |kg CH2O ha-1|
    MREST         Total gross maintenance respiration                N    |kg CH2O ha-1|
    CTRAT         Total crop transpiration                           N    cm
    HI            Harvest Index (only calculated during              N    -
                  `finalize()`)
    DOF           Date representing the day of finish of the crop    N    -
                  simulation.
    FINISH_TYPE   String representing the reason for finishing the   N    -
                  simulation: maturity, harvest, leave death, etc.
    ============  ================================================= ==== ===============

 
     **Rate variables:**

    =======  ================================================ ==== =============
     Name     Description                                      Pbl      Unit
    =======  ================================================ ==== =============
    GASS     Assimilation rate corrected for water stress       N  |kg CH2O ha-1 d-1|
    PGASS    Potential assimilation rate                        N  |kg CH2O ha-1 d-1|
    MRES     Actual maintenance respiration rate, taking into
             account that MRES <= GASS.                         N  |kg CH2O ha-1 d-1|
    PMRES    Potential maintenance respiration rate             N  |kg CH2O ha-1 d-1|
    ASRC     Net available assimilates (GASS - MRES)            N  |kg CH2O ha-1 d-1|
    DMI      Total dry matter increase, calculated as ASRC
             times a weighted conversion efficieny.             Y  |kg ha-1 d-1|
    ADMI     Aboveground dry matter increase                    Y  |kg ha-1 d-1|
    =======  ================================================ ==== =============

    """
    
    # sub-model components for crop simulation
    pheno = Instance(SimulationObject)
    part = Instance(SimulationObject)
    assim = Instance(SimulationObject)
    mres = Instance(SimulationObject)
    evtra = Instance(SimulationObject)
    lv_dynamics = Instance(SimulationObject)
    st_dynamics = Instance(SimulationObject)
    ro_dynamics = Instance(SimulationObject)
    so_dynamics = Instance(SimulationObject)
    npk_crop_dynamics = Instance(SimulationObject)
    npk_soil_dynamics = Instance(SimulationObject)
    npk_stress = Instance(SimulationObject)
        
    # Parameters, rates and states which are relevant at the main crop
    # simulation level
    class Parameters(ParamTemplate):
        CVL = Float(-99.)
        CVO = Float(-99.)
        CVR = Float(-99.)
        CVS = Float(-99.)

    class StateVariables(StatesTemplate):
        TAGP = Float(-99.)
        GASST = Float(-99.)
        MREST = Float(-99.)
        CTRAT = Float(-99.) # Crop total transpiration
        HI = Float(-99.)
        DOF = Instance(datetime.date)
        FINISH_TYPE = Instance(str)

    class RateVariables(RatesTemplate):
        GASS = Float(-99.)
        PGASS = Float(-99.)
        MRES = Float(-99.)
        PMRES = Float(-99.)
        ASRC = Float(-99.)
        DMI = Float(-99.)
        ADMI = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE model instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["DMI", "ADMI"])
        self.kiosk = kiosk
        
        # Initialize components of the crop
        self.pheno = Phenology(day, kiosk,  parvalues)
        self.part = Partitioning(day, kiosk, parvalues)
        self.assim = Assimilation(day, kiosk, parvalues)
        self.mres = MaintenanceRespiration(day, kiosk, parvalues)
        self.evtra = Evapotranspiration(day, kiosk, parvalues)
        self.ro_dynamics = Root_Dynamics(day, kiosk, parvalues)
        self.st_dynamics = Stem_Dynamics(day, kiosk, parvalues)
        self.so_dynamics = Storage_Organ_Dynamics(day, kiosk, parvalues)
        self.lv_dynamics = Leaf_Dynamics(day, kiosk, parvalues)
        # Added for book keeping of N/P/K in crop and soil
        self.npk_crop_dynamics = NPK_crop(day, kiosk, parvalues)
        self.npk_soil_dynamics = NPK_soil(day, kiosk, parvalues)
        self.npk_stress = NPK_Stress(day, kiosk, parvalues)
        

        # Initial total (living+dead) above-ground biomass of the crop
        TAGP = self.kiosk["TWLV"] + \
               self.kiosk["TWST"] + \
               self.kiosk["TWSO"]
        
        self.states = self.StateVariables(kiosk,
                                          publish=["TAGP","GASST","MREST","HI"],
                                          TAGP=TAGP, GASST=0.0, MREST=0.0,
                                          CTRAT=0.0, HI=0.0,
                                          DOF=None, FINISH_TYPE=None)

        # Check partitioning of TDWI over plant organs
        checksum = parvalues["TDWI"] - self.states.TAGP - self.kiosk["TWRT"]
        if abs(checksum) > 0.0001:
            msg = "Error in partitioning of initial biomass (TDWI)!"
            raise exc.PartitioningError(msg)
            
        # assign handler for CROP_FINISH signal
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)

    @staticmethod
    def _check_carbon_balance(day, DMI, GASS, MRES, CVF, pf):
        (FR, FL, FS, FO) = pf
        checksum = (GASS - MRES - (FR+(FL+FS+FO)*(1.-FR)) * DMI/CVF) * \
                    1./(max(0.0001,GASS))
        if abs(checksum) >= 0.0001:
            msg = "Carbon flows not balanced on day %s\n" % day
            msg += "Checksum: %f, GASS: %f, MRES: %f\n"  (checksum, GASS, MRES)
            msg += "FR,L,S,O: %5.3f,%5.3f,%5.3f,%5.3f, DMI: %f, CVF: %f\n" % \
                   (FR,FL,FS,FO,DMI,CVF)
            raise exc.CarbonBalanceError(msg)

    @prepare_rates
    def calc_rates(self, day, drv):
        params = self.params
        rates  = self.rates

        # Phenology
        self.pheno.calc_rates(day, drv)
        crop_stage = self.pheno.get_variable("STAGE")

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if crop_stage == "emerging":
            return

        # Potential assimilation
        rates.PGASS = self.assim(day, drv)
        
        # (evapo)transpiration rates
        self.evtra(day, drv)

        # nutrient status and reduction factor
        NNI, NPKI, NPKREF = self.npk_stress(day, drv)

        # water stress reduction
        TRA = self.kiosk["TRA"]
        TRAMX = self.kiosk["TRAMX"]
        
        # Select minimum of nutrient and water stress
        reduction = min(NPKREF, TRA/TRAMX)

        rates.GASS = rates.PGASS * reduction

        # Respiration
        rates.PMRES = self.mres(day, drv)
        rates.MRES = min(rates.GASS, rates.PMRES)

        # Net available assimilates
        rates.ASRC = rates.GASS - rates.MRES

        # DM partitioning factors (pf), conversion factor (CVF),
        # dry matter increase (DMI) and check on carbon balance
        pf = self.part.calc_rates(day, drv)
        CVF = 1./((pf.FL/params.CVL + pf.FS/params.CVS + pf.FO/params.CVO) *
                  (1.-pf.FR) + pf.FR/params.CVR)
        rates.DMI = CVF * rates.ASRC
        self._check_carbon_balance(day, rates.DMI, rates.GASS, rates.MRES,
                                   CVF, pf)

        # distribution over plant organ

        # Below-ground dry matter increase and root dynamics
        self.ro_dynamics.calc_rates(day, drv)
        # Aboveground dry matter increase and distribution over stems,
        # leaves, organs
        rates.ADMI = (1. - pf.FR) * rates.DMI
        self.st_dynamics.calc_rates(day, drv)
        self.so_dynamics.calc_rates(day, drv)
        self.lv_dynamics.calc_rates(day, drv)
        
        # Update nutrient rates in crop and soil
        self.npk_crop_dynamics.calc_rates(day, drv)
        self.npk_soil_dynamics.calc_rates(day, drv)
        
    @prepare_states
    def integrate(self, day, delt=1.0):
        rates = self.rates
        states = self.states

        # crop stage before integration
        crop_stage = self.pheno.get_variable("STAGE")

        # Phenology
        self.pheno.integrate(day, delt)

        # if before emergence there is no need to continue
        # because only the phenology is running.
        # Just run a touch() to to ensure that all state variables are available
        # in the kiosk
        if crop_stage == "emerging":
            self.touch()
            return

        # Partitioning
        self.part.integrate(day, delt)
        
        # Integrate states on leaves, storage organs, stems and roots
        self.ro_dynamics.integrate(day, delt)
        self.so_dynamics.integrate(day, delt)
        self.st_dynamics.integrate(day, delt)
        self.lv_dynamics.integrate(day, delt)

        # Update nutrient states in crop and soil
        self.npk_crop_dynamics.integrate(day, delt)
        self.npk_soil_dynamics.integrate(day, delt)
        

        # Integrate total (living+dead) above-ground biomass of the crop
        states.TAGP = self.kiosk["TWLV"] + \
                      self.kiosk["TWST"] + \
                      self.kiosk["TWSO"]

        # total gross assimilation and maintenance respiration 
        states.GASST += rates.GASS
        states.MREST += rates.MRES
        
        # total crop transpiration (CTRAT)
        states.CTRAT += self.kiosk["TRA"]
        
    @prepare_states
    def finalize(self, day):

        # Calculate Harvest Index
        if self.states.TAGP > 0:
            self.states.HI = self.kiosk["TWSO"]/self.states.TAGP
        else:
            msg = "Cannot calculate Harvest Index because TAGP=0"
            self.logger.warning(msg)
            self.states.HI = -1.
        
        SimulationObject.finalize(self, day)

    def _on_CROP_FINISH(self, day, finish_type=None):
        """Handler for setting day of finish (DOF) and reason for
        crop finishing (FINISH).
        """
        self._for_finalize["DOF"] = day
        self._for_finalize["FINISH_TYPE"] = finish_type
