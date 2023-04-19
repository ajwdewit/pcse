# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from math import exp
from collections import deque
from array import array
import numpy as np

from ..traitlets import Float, Int, Instance
from ..decorators import prepare_rates, prepare_states
from ..util import limit, AfgenTrait
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, \
     SimulationObject
from .. import signals


class WOFOST72_Leaf_Dynamics(SimulationObject):
    """Leaf dynamics for the WOFOST crop model.
    
    Implementation of biomass partitioning to leaves, growth and senenscence
    of leaves. WOFOST keeps track of the biomass that has been partitioned to
    the leaves for each day (variable `LV`), which is called a leaf class).
    For each leaf class the leaf age (variable 'LVAGE') and specific leaf area
    (variable `SLA`) are also registered. Total living leaf biomass is
    calculated by summing the biomass values for all leaf classes. Similarly,
    leaf area is calculated by summing leaf biomass times specific leaf area
    (`LV` * `SLA`).
    
    Senescense of the leaves can occur as a result of physiological age,
    drought stress or self-shading.
       
    *Simulation parameters* (provide in cropdata dictionary)
    
    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    RGRLAI   Maximum relative increase in LAI.              SCr     ha ha-1 d-1
    SPAN     Life span of leaves growing at 35 Celsius      SCr     |d|
    TBASE    Lower threshold temp. for ageing of leaves     SCr     |C|
    PERDL    Max. relative death rate of leaves due to      SCr
             water stress
    TDWI     Initial total crop dry weight                  SCr     |kg ha-1|
    KDIFTB   Extinction coefficient for diffuse visible     TCr
             light as function of DVS
    SLATB    Specific leaf area as a function of DVS        TCr     |ha kg-1|
    =======  ============================================= =======  ============

    *State variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    LV       Leaf biomass per leaf class                        N    |kg ha-1|
    SLA      Specific leaf area per leaf class                  N    |ha kg-1|
    LVAGE    Leaf age per leaf class                            N    |d|
    LVSUM    Sum of LV                                          N    |kg ha-1|
    LAIEM    LAI at emergence                                   N    -
    LASUM    Total leaf area as sum of LV*SLA,                  N    -
             not including stem and pod area                    N
    LAIEXP   LAI value under theoretical exponential growth     N    -
    LAIMAX   Maximum LAI reached during growth cycle            N    -
    LAI      Leaf area index, including stem and pod area       Y    -
    WLV      Dry weight of living leaves                        Y    |kg ha-1|
    DWLV     Dry weight of dead leaves                          N    |kg ha-1|
    TWLV     Dry weight of total leaves (living + dead)         Y    |kg ha-1|
    =======  ================================================= ==== ============


    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    GRLV     Growth rate leaves                                 N   |kg ha-1 d-1|
    DSLV1    Death rate leaves due to water stress              N   |kg ha-1 d-1|
    DSLV2    Death rate leaves due to self-shading              N   |kg ha-1 d-1|
    DSLV3    Death rate leaves due to frost kill                N   |kg ha-1 d-1|
    DSLV     Maximum of DLSV1, DSLV2, DSLV3                     N   |kg ha-1 d-1|
    DALV     Death rate leaves due to aging.                    N   |kg ha-1 d-1|
    DRLV     Death rate leaves as a combination of DSLV and     N   |kg ha-1 d-1|
             DALV
    SLAT     Specific leaf area for current time step,          N   |ha kg-1|
             adjusted for source/sink limited leaf expansion
             rate.
    FYSAGE   Increase in physiological leaf age                 N   -
    GLAIEX   Sink-limited leaf expansion rate (exponential      N   |ha ha-1 d-1|
             curve) 
    GLASOL   Source-limited leaf expansion rate (biomass        N   |ha ha-1 d-1|
             increase)
    =======  ================================================= ==== ============

    
    *External dependencies:*
    
    ======== ============================== =============================== ===========
     Name     Description                         Provided by               Unit
    ======== ============================== =============================== ===========
    DVS      Crop development stage         DVS_Phenology                    - 
    FL       Fraction biomass to leaves     DVS_Partitioning                 -
    FR       Fraction biomass to roots      DVS_Partitioning                 -
    SAI      Stem area index                WOFOST_Stem_Dynamics             -
    PAI      Pod area index                 WOFOST_Storage_Organ_Dynamics    -
    TRA      Transpiration rate             Evapotranspiration              |cm day-1|
    TRAMX    Maximum transpiration rate     Evapotranspiration              |cm day-1| 
    ADMI     Above-ground dry matter        CropSimulation                  |kg ha-1 d-1|
             increase
    RF_FROST Reduction factor frost kill    FROSTOL                          -
    ======== ============================== =============================== ===========
    """

    class Parameters(ParamTemplate):
        RGRLAI = Float(-99.)
        SPAN   = Float(-99.)
        TBASE  = Float(-99.)
        PERDL  = Float(-99.)
        TDWI   = Float(-99.)
        SLATB  = AfgenTrait()
        KDIFTB = AfgenTrait()

    class StateVariables(StatesTemplate):
        LV     = Instance(deque)
        SLA    = Instance(deque)
        LVAGE  = Instance(deque)
        LAIEM  = Float(-99.)
        LASUM  = Float(-99.)
        LAIEXP = Float(-99.)
        LAIMAX = Float(-99.)
        LAI    = Float(-99.)
        WLV    = Float(-99.)
        DWLV   = Float(-99.)
        TWLV   = Float(-99.)

    class RateVariables(RatesTemplate):
        GRLV  = Float(-99.)
        DSLV1 = Float(-99.)
        DSLV2 = Float(-99.)
        DSLV3 = Float(-99.)
        DSLV  = Float(-99.)
        DALV  = Float(-99.)
        DRLV  = Float(-99.)
        SLAT  = Float(-99.)
        FYSAGE = Float(-99.)
        GLAIEX = Float(-99.)
        GLASOL = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)
        self.rates  = self.RateVariables(kiosk)

        p = self.params
        k = self.kiosk

        # Initial leaf biomass
        WLV  = (p.TDWI * (1-k.FR)) * k.FL

        # Initial values for leaf area
        LAIEM  = WLV * p.SLATB(k.DVS)

        s = dict(
            WLV=WLV,
            DWLV = 0.,
            TWLV = WLV,
            # First leaf class (SLA, age and weight)
            SLA = deque([p.SLATB(k.DVS)]),
            LVAGE = deque([0.]),
            LV = deque([WLV]),
            LAIEM = LAIEM,
            LASUM = LAIEM,
            LAIEXP = LAIEM,
            LAIMAX = LAIEM,
            LAI = LAIEM + k.SAI + k.PAI
        )

        self.states = self.StateVariables(kiosk, publish=["LAI","TWLV","WLV"], **s)

    def _calc_LAI(self):
        # Total leaf area Index as sum of leaf, pod and stem area
        SAI = self.kiosk["SAI"]
        PAI = self.kiosk["PAI"]
        return self.states.LASUM + SAI + PAI

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        # Growth rate leaves
        # weight of new leaves
        r.GRLV = k.ADMI * k.FL

        # death of leaves due to water/oxygen stress
        r.DSLV1 = s.WLV * (1. - k.RFTRA) * p.PERDL

        # death due to self shading cause by high LAI
        LAICR = 3.2/p.KDIFTB(k.DVS)
        r.DSLV2 = s.WLV * limit(0., 0.03, 0.03*(s.LAI-LAICR)/LAICR)

        # Death of leaves due to frost damage as determined by
        # Reduction Factor Frost "RF_FROST"
        if "RF_FROST" in self.kiosk:
            r.DSLV3 = s.WLV * k.RF_FROST
        else:
            r.DSLV3 = 0.

        # leaf death equals maximum of water stress, shading and frost
        r.DSLV = max(r.DSLV1, r.DSLV2, r.DSLV3)

        # Determine how much leaf biomass classes have to die in states.LV,
        # given the life span > SPAN, these classes will be accumulated
        # in DALV.
        # Note that the actual leaf death is imposed on the array LV during the
        # state integration step.
        DALV = 0.0
        for lv, lvage in zip(s.LV, s.LVAGE):
            if lvage > p.SPAN:
                DALV += lv
        r.DALV = DALV

        # Total death rate leaves
        r.DRLV = max(r.DSLV, r.DALV)

        # physiologic ageing of leaves per time step
        r.FYSAGE = max(0., (drv.TEMP - p.TBASE)/(35. - p.TBASE))

        # specific leaf area of leaves per time step
        r.SLAT = p.SLATB(k.DVS)

        # leaf area not to exceed exponential growth curve
        if s.LAIEXP < 6.:
            DTEFF = max(0., drv.TEMP-p.TBASE)
            r.GLAIEX = s.LAIEXP * p.RGRLAI * DTEFF
            # source-limited increase in leaf area
            r.GLASOL = r.GRLV * r.SLAT
            # sink-limited increase in leaf area
            GLA = min(r.GLAIEX, r.GLASOL)
            # adjustment of specific leaf area of the youngest leaf class
            if r.GRLV > 0.:
                r.SLAT = GLA/r.GRLV

    @prepare_states
    def integrate(self, day, delt=1.0):
        p = self.params
        r = self.rates
        s = self.states

        # --------- leave death ---------
        tLV = array('d', s.LV)
        tSLA = array('d', s.SLA)
        tLVAGE = array('d', s.LVAGE)
        tDRLV = r.DRLV

        # leaf death is imposed on leaves by removing leave classes from the
        # right side of the deque. 
        for LVweight in reversed(s.LV):
            if tDRLV > 0.:
                if tDRLV >= LVweight: # remove complete leaf class from deque
                    tDRLV -= LVweight
                    tLV.pop()
                    tLVAGE.pop()
                    tSLA.pop()
                else: # Decrease value of oldest (rightmost) leave class
                    tLV[-1] -= tDRLV
                    tDRLV = 0.
            else:
                break

        # Integration of physiological age
        tLVAGE = deque([age + r.FYSAGE for age in tLVAGE])
        tLV = deque(tLV)
        tSLA = deque(tSLA)

        # --------- leave growth ---------
        # new leaves in class 1
        tLV.appendleft(r.GRLV)
        tSLA.appendleft(r.SLAT)
        tLVAGE.appendleft(0.)

        # calculation of new leaf area
        s.LASUM = sum([lv*sla for lv, sla in zip(tLV, tSLA)])
        s.LAI = self._calc_LAI()
        s.LAIMAX = max(s.LAI, s.LAIMAX)

        # exponential growth curve
        s.LAIEXP += r.GLAIEX

        # Update leaf biomass states
        s.WLV  = sum(tLV)
        s.DWLV += r.DRLV
        s.TWLV = s.WLV + s.DWLV

        # Store final leaf biomass deques
        self.states.LV = tLV
        self.states.SLA = tSLA
        self.states.LVAGE = tLVAGE

    @prepare_states
    def _set_variable_LAI(self, nLAI):
        """Updates the value of LAI to the new value provided as input.

        Related state variables will be updated as well and the increments
        to all adjusted state variables will be returned as a dict.
        """
        s = self.states
        k = self.kiosk

        # Store old values of states
        oWLV = s.WLV
        oLAI = s.LAI
        oTWLV = s.TWLV
        oLASUM = s.LASUM

        # Reduce oLAI for pod and stem area. SAI and PAI will not be adjusted
        # because this is often only a small component of the total leaf
        # area. For all current crop files in WOFOST SPA and SSA are zero
        # anyway
        SAI = k.SAI
        PAI = k.PAI
        adj_nLAI = max(nLAI - SAI - PAI, 0.)
        adj_oLAI = max(oLAI - SAI - PAI, 0.)

        # LAI Adjustment factor for leaf biomass LV (rLAI)
        if adj_oLAI > 0:
            rLAI = adj_nLAI/adj_oLAI
            LV = [lv*rLAI for lv in s.LV]
        # If adj_oLAI == 0 then add the leave biomass directly to the
        # youngest leave age class (LV[0])
        else:
            LV = [nLAI/s.SLA[0]]

        s.LASUM = sum([lv*sla for lv, sla in zip(LV, s.SLA)])
        s.LV = deque(LV)
        s.LAI = self._calc_LAI()
        s.WLV = sum(s.LV)
        s.TWLV = s.WLV + s.DWLV

        increments = {"LAI": s.LAI - oLAI,
                      "LAISUM":s.LASUM - oLASUM,
                      "WLV": s.WLV - oWLV,
                      "TWLV": s.TWLV - oTWLV}
        return increments


class CSDM_Leaf_Dynamics(SimulationObject):
    """Leaf dynamics according to the Canopy Structure Dynamic Model.
    
    The only difference is that in the real CSDM the temperature sum is the
    driving variable, while in this case it is simply the day number since \
    the start of the model.
    
    Reference:
    Koetz et al. 2005. Use of coupled canopy structure dynamic and radiative
    transfer models to estimate biophysical canopy characteristics.
    Remote Sensing of Environment. Volume 95, Issue 1, 15 March 2005,
    Pages 115-124. http://dx.doi.org/10.1016/j.rse.2004.11.017

        
    For plotting the CSDM model with GNUPLOT the following example code can be used:

        td = 150
        CSDM_MAX = 5.
        CSDM_MIN = 0.15
        CSDM_A = 0.085
        CSDM_B = 0.045
        CSDM_T1 = int(td/3.)
        CSDM_T2 = td

        set xrange [0:200]
        set yrange [-1:8]
        plot CSDM_MIN + CSDM_MAX*(1./(1. + exp(-CSDM_B*(x - CSDM_T1)))**2 - exp(CSDM_A*(x - CSDM_T2)))
    
    """

    class Parameters(ParamTemplate):
        CSDM_MAX = Float()
        CSDM_MIN = Float()
        CSDM_A = Float()
        CSDM_B = Float()
        CSDM_T1 = Float()
        CSDM_T2 = Float()

    class StateVariable(StatesTemplate):
        LAI = Float()
        DAYNR = Int()
        LAIMAX = Float()

    def _CSDM(self, daynr):
        """Returns the LAI value depending on the day number (daynr)
        and the parameters of the CSDM model.
        """
        p = self.params

        LAI_growth = 1./(1. + exp(-p.CSDM_B*(daynr - p.CSDM_T1)))**2
        LAI_senescence = -exp(p.CSDM_A*(daynr - p.CSDM_T2))
        LAI = p.CSDM_MIN + p.CSDM_MAX*(LAI_growth + LAI_senescence)

        # do not allow LAI lower then CSDM_MIN
        if LAI < p.CSDM_MIN:
            msg = ("LAI of CSDM model smaller then lower LAI limit "+
                   "(CSDM_MIN)! Adjusting LAI to CSDM_MIN.")
            self.logger.warn(msg)
            LAI = max(p.CSDM_MIN, LAI)

        return LAI

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.params = self.Parameters(parvalues)

        # calculate LAI on day 1 from CSDM
        LAI = self._CSDM(1)
        self.states = self.StateVariable(kiosk, LAI=LAI, DAYNR=1,
                                         LAIMAX=self.params.CSDM_MIN,
                                         publish="LAI")

    @prepare_rates
    def calc_rates(self, day, drv):
        pass

    @prepare_states
    def integrate(self, day, delt=1.0):

        self.states.DAYNR += 1
        self.states.LAI = self._CSDM(self.states.DAYNR)
        if self.states.LAI > self.states.LAIMAX:
            self.states.LAIMAX = self.states.LAI

        if self.states.DAYNR > self.params.CSDM_T2:
            self._send_signal(signal=signals.crop_finish, day=day,
                              finish_type="Canopy died according to CSDM leaf model.",
                              crop_delete=True)


class WOFOST81_Leaf_Dynamics(SimulationObject):
    """Leaf dynamics for the WOFOST crop model including leaf response to
    N stress.

    # TODO HB 20220405: This function was changed quite a bit and needs redocumentation.

    Implementation of biomass partitioning to leaves, growth and senenscence
    of leaves. WOFOST keeps track of the biomass that has been partitioned to
    the leaves for each day (variable `LV`), which is called a leaf class).
    For each leaf class the leaf age (variable 'LVAGE') and specific leaf area
    are (variable `SLA`) are also registered. Total living leaf biomass
    is calculated by summing the biomass values for all leaf classes. Similarly,
    leaf area is calculated by summing leaf biomass times specific leaf area
    (`LV` * `SLA`).

    Senescense of the leaves can occur as a result of physiological age,
    drought stress, nutrient stress or self-shading.

    Finally, leaf expansion (SLA) can be influenced by nutrient stress.

    *Simulation parameters*

    =======  ============================================= =======  ============
     Name     Description                                   Type     Unit
    =======  ============================================= =======  ============
    RGRLAI   Maximum relative increase in LAI.              SCr     ha ha-1 d-1
    SPAN     Life span of leaves growing at 35 Celsius      SCr     |d|
    TBASE    Lower threshold temp. for ageing of leaves     SCr     |C|
    PERDL    Max. relative death rate of leaves due to      SCr
             water stress
    TDWI     Initial total crop dry weight                  SCr     |kg ha-1|
    KDIFTB   Extinction coefficient for diffuse visible     TCr
             light as function of DVS
    SLATB    Specific leaf area as a function of DVS        TCr     |ha kg-1|
    RDRNS    max. relative death rate of leaves due to      TCr         -
             nutrient N stress
    NLAI     coefficient for the reduction due to           TCr         -
             nutrient N stress of the LAI increase
             (during juvenile phase).
    NSLA     Coefficient for the effect of nutrient NPK     TCr         -
             stress on SLA reduction
    RDRNS    Max. relative death rate of leaves due to      TCr         -
             nutrient N stress
    =======  ============================================= =======  ============

    *State variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    LV       Leaf biomass per leaf class                        N    |kg ha-1|
    SLA      Specific leaf area per leaf class                  N    |ha kg-1|
    LVAGE    Leaf age per leaf class                            N    |d|
    LVSUM    Sum of LV                                          N    |kg ha-1|
    LAIEM    LAI at emergence                                   N    -
    LASUM    Total leaf area as sum of LV*SLA,                  N    -
             not including stem and pod area                    N
    LAIEXP   LAI value under theoretical exponential growth     N    -
    LAIMAX   Maximum LAI reached during growth cycle            N    -
    LAI      Leaf area index, including stem and pod area       Y    -
    WLV      Dry weight of living leaves                        Y    |kg ha-1|
    DWLV     Dry weight of dead leaves                          N    |kg ha-1|
    TWLV     Dry weight of total leaves (living + dead)         Y    |kg ha-1|
    =======  ================================================= ==== ============


    *Rate variables*

    =======  ================================================= ==== ============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ============
    GRLV     Growth rate leaves                                 N   |kg ha-1 d-1|
    DSLV1    Death rate leaves due to water stress              N   |kg ha-1 d-1|
    DSLV2    Death rate leaves due to self-shading              N   |kg ha-1 d-1|
    DSLV3    Death rate leaves due to frost kill                N   |kg ha-1 d-1|
    DSLV4    Death rate leaves due to nutrient stress           N   |kg ha-1 d-1|
    DSLV     Maximum of DLSV1, DSLV2, DSLV3                     N   |kg ha-1 d-1|
    DALV     Death rate leaves due to aging.                    N   |kg ha-1 d-1|
    DRLV     Death rate leaves as a combination of DSLV and     N   |kg ha-1 d-1|
             DALV
    SLAT     Specific leaf area for current time step,          N   |ha kg-1|
             adjusted for source/sink limited leaf expansion
             rate.
    FYSAGE   Increase in physiological leaf age                 N   -
    GLAIEX   Sink-limited leaf expansion rate (exponential      N   |ha ha-1 d-1|
             curve)
    GLASOL   Source-limited leaf expansion rate (biomass        N   |ha ha-1 d-1|
             increase)
    =======  ================================================= ==== ============


    *External dependencies:*

    ======== ============================== =============================== ===========
     Name     Description                         Provided by               Unit
    ======== ============================== =============================== ===========
    DVS      Crop development stage         DVS_Phenology                    -
    FL       Fraction biomass to leaves     DVS_Partitioning                 -
    FR       Fraction biomass to roots      DVS_Partitioning                 -
    SAI      Stem area index                WOFOST_Stem_Dynamics             -
    PAI      Pod area index                 WOFOST_Storage_Organ_Dynamics    -
    TRA      Transpiration rate             Evapotranspiration              |cm day-1|
    TRAMX    Maximum transpiration rate     Evapotranspiration              |cm day-1|
    ADMI     Above-ground dry matter        CropSimulation                  |kg ha-1 d-1|
             increase
    RF_FROST Reduction factor frost kill    FROSTOL                          -
    ======== ============================== =============================== ===========
    """

    class Parameters(ParamTemplate):
        RGRLAI = Float(-99.)
        SPAN = Float(-99.)
        TBASE = Float(-99.)
        PERDL = Float(-99.)
        TDWI = Float(-99.)
        SLATB = AfgenTrait()
        KDIFTB = AfgenTrait()

    class StateVariables(StatesTemplate):
        LV = Instance(deque)
        SLA = Instance(deque)
        LVAGE = Instance(deque)
        LAIEM = Float(-99.)
        LASUM = Float(-99.)
        LAIEXP = Float(-99.)
        LAIMAX = Float(-99.)
        LAI = Float(-99.)
        WLV = Float(-99.)
        DWLV = Float(-99.)
        TWLV = Float(-99.)

    class RateVariables(RatesTemplate):
        GRLV = Float(-99.)
        DSLV1 = Float(-99.)
        DSLV2 = Float(-99.)
        DSLV3 = Float(-99.)
        DSLV4 = Float(-99.)
        DSLV = Float(-99.)
        DALV = Float(-99.)
        DRLV = Float(-99.)
        SLAT = Float(-99.)
        FYSAGE = Float(-99.)
        GLAIEX = Float(-99.)
        GLASOL = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with WOFOST cropdata key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["DRLV", "GRLV"])

        # CALCULATE INITIAL STATE VARIABLES
        p = self.params
        k = self.kiosk

        # Initial leaf biomass
        WLV = (p.TDWI * (1-k.FR)) * k.FL

        # # Initial values for leaf area
        LAIEM = WLV * p.SLATB(k.DVS)

        s = dict(
            WLV=WLV,
            DWLV = 0.,
            TWLV = WLV,
            # First leaf class (SLA, age and weight)
            SLA = deque([p.SLATB(k.DVS)]),
            LVAGE = deque([0.]),
            LV = deque([WLV]),
            LAIEM = LAIEM,
            LASUM = LAIEM,
            LAIEXP = LAIEM,
            LAIMAX = LAIEM,
            LAI = LAIEM + k.SAI + k.PAI
        )

        # Initialize StateVariables object
        self.states = self.StateVariables(kiosk, publish=["LAI", "TWLV", "WLV"], **s)

    def _calc_LAI(self):
        # Total leaf area Index as sum of leaf, pod and stem area
        k = self.kiosk
        return self.states.LASUM + k.SAI + k.PAI

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        # Growth rate leaves
        # weight of new leaves
        r.GRLV = k.ADMI * k.FL

        # death of leaves due to water/oxygen stress
        r.DSLV1 = s.WLV * (1.-k.RFTRA) * p.PERDL

        # death due to self shading cause by high LAI
        LAICR = 3.2/p.KDIFTB(k.DVS)
        r.DSLV2 = s.WLV * limit(0., 0.03, 0.03*(s.LAI-LAICR)/LAICR)

        # Death of leaves due to frost damage as determined by
        # Reduction Factor Frost "RF_FROST"
        if "RF_FROST" in k:
            r.DSLV3 = s.WLV * k.RF_FROST
        else:
            r.DSLV3 = 0.

        # Determine how much leaf biomass classes have to die in states.LV,
        # given the life span > SPAN, these classes will be accumulated
        # in DALV.
        # Note that the actual leaf death is imposed on the array LV during the
        # state integration step.
        DALV = 0.0
        for lv, lvage in zip(s.LV, s.LVAGE):
            if lvage > p.SPAN:
                DALV += lv
        r.DALV = DALV

        # Suggestion Allard accelerate ageing:
        r.DSLV = max(r.DSLV1, r.DSLV2, r.DSLV3)
        r.DALV = min(DALV * k.NSLLV, k.WLV)

        # Total death rate leaves
        r.DRLV = max(r.DSLV, r.DALV)

        # physiologic ageing of leaves per time step
        r.FYSAGE = max(0., (drv.TEMP - p.TBASE)/(35. - p.TBASE))

        ## specific leaf area of leaves per time step
        r.SLAT = p.SLATB(k.DVS)

        # leaf area not to exceed exponential growth curve
        if s.LAIEXP < 6.:
            DTEFF = max(0., drv.TEMP-p.TBASE)

            # added IS
            # Nutrient and water stress during juvenile stage:
            if k.DVS < 0.2 and s.LAI < 0.75:
                factor = k.RFTRA * k.RFRGRL
            else:
                factor = 1.

            r.GLAIEX = s.LAIEXP * p.RGRLAI * DTEFF * factor
            # source-limited increase in leaf area
            r.GLASOL = r.GRLV * r.SLAT
            # sink-limited increase in leaf area
            GLA = min(r.GLAIEX, r.GLASOL)
            # adjustment of specific leaf area of the youngest leaf class
            if r.GRLV > 0.:
                r.SLAT = GLA/r.GRLV

    @prepare_states
    def integrate(self, day, delt=1.0):
        p = self.params
        r = self.rates
        s = self.states
        k = self.kiosk

        # --------- leave death ---------
        tLV = array('d', s.LV)
        tSLA = array('d', s.SLA)
        tLVAGE = array('d', s.LVAGE)
        tDRLV = r.DRLV

        # leaf death is imposed on leaves by removing leave classes from the
        # right side of the deque.
        for LVweight in reversed(s.LV):
            if tDRLV > 0.:
                if tDRLV >= LVweight: # remove complete leaf class from deque
                    tDRLV -= LVweight
                    tLV.pop()
                    tLVAGE.pop()
                    tSLA.pop()
                else: # Decrease value of oldest (rightmost) leave class
                    tLV[-1] -= tDRLV
                    tDRLV = 0.
            else:
                break

        # Integration of physiological age
        tLVAGE = deque([age + r.FYSAGE for age in tLVAGE])

        # Uniformly reduce leaf biomass in case of reallocation
        if k.REALLOC_LV > 0:
            sumLV = sum(tLV)
            if k.REALLOC_LV < sumLV:
                ReductionFactorLV = (sumLV - k.REALLOC_LV)/sumLV
                tLV = np.array(tLV) * ReductionFactorLV

        tLV = deque(tLV)
        tSLA = deque(tSLA)

        # --------- leave growth ---------
        # new leaves in class 1
        tLV.appendleft(r.GRLV)
        tSLA.appendleft(r.SLAT)
        tLVAGE.appendleft(0.)

        # calculation of new leaf area
        s.LASUM = sum([lv*sla for lv, sla in zip(tLV, tSLA)])
        s.LAI = self._calc_LAI()
        s.LAIMAX = max(s.LAI, s.LAIMAX)

        # exponential growth curve
        s.LAIEXP += r.GLAIEX

        # Update leaf biomass states
        s.WLV  = sum(tLV)
        s.DWLV += r.DRLV
        s.TWLV = s.WLV + s.DWLV

        # Store final leaf biomass deques
        self.states.LV = tLV
        self.states.SLA = tSLA
        self.states.LVAGE = tLVAGE
