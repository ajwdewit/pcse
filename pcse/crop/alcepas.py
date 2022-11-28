# -*- coding: utf-8 -*-
# Copyright (c) 2004-2017 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), October 2017
"""Basic routines for ALCEPAS onion model
"""
from __future__ import print_function
from math import exp
from collections import deque
from array import array
import datetime as dt

from ..traitlets import Instance, Float, Enum, Unicode
from .assimilation import totass
from ..util import limit, astro, doy, daylength, AfgenTrait
from ..base import ParamTemplate, StatesTemplate, RatesTemplate, SimulationObject
from ..decorators import prepare_rates, prepare_states
from .. import signals
from .. import exceptions as exc


class Respiration(SimulationObject):

    class Parameters(ParamTemplate):
        Q10 = Float()
        MSOTB = AfgenTrait()
        MLVTB = AfgenTrait()
        MRTTB = AfgenTrait()

    class RateVariables(RatesTemplate):
        MAINT = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="MAINT")

    @prepare_rates
    def __call__(self, day, drv):
        r = self.rates
        k = self.kiosk
        p = self.params

        MAINSO = p.MSOTB(k.DVS)
        MAINLV = p.MLVTB(k.DVS)
        MAINRT = p.MRTTB(k.DVS)
        MAINTS = MAINLV * k.WLV + MAINRT * k.WRT + MAINSO * k.WSO

        TEFF = p.Q10 ** ((drv.TEMP - 20.) / 10.)
        MNDVS = 1.0
        r.MAINT = min(k.GPHOT, MAINTS * TEFF * MNDVS)
        return r.MAINT


class Assimilation(SimulationObject):

    class Parameters(ParamTemplate):
        AMX = Float()
        EFF = Float()
        KDIF = Float()
        AMDVST = AfgenTrait()
        AMTMPT = AfgenTrait()

    class RateVariables(RatesTemplate):
        GPHOT = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish="GPHOT")

    @prepare_rates
    def __call__(self, day, drv):
        r = self.rates
        k = self.kiosk
        p = self.params

        a = astro(day, drv.LAT, drv.IRRAD)
        AMDVS = p.AMDVST(k.DVS)
        AMTMP = p.AMTMPT(drv.DTEMP)
        AMAX = p.AMX * AMDVS * AMTMP
        DTGA = totass(a.DAYL, AMAX, p.EFF, k.LAI, p.KDIF, drv.IRRAD, a.DIFPP, a.DSINBE, a.SINLD, a.COSLD)
        r.GPHOT = DTGA * 30./44.
        return r.GPHOT


class Partitioning(SimulationObject):

    class Parameters(ParamTemplate):
        FLVTB = AfgenTrait()
        FSHTB = AfgenTrait()

    class RateVariables(RatesTemplate):
        FSH = Float()
        FLV = Float()
        FSO = Float()
        FRT = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["FSH","FLV","FRT","FSO"])

    @prepare_rates
    def __call__(self, day, drv):
        r = self.rates
        p = self.params
        k = self.kiosk
        r.FSH = p.FSHTB(k.DVS)
        r.FRT = 1. - r.FSH
        r.FLV = p.FLVTB(k.DVS)
        r.FSO = 1. - r.FLV


class Phenology(SimulationObject):

    class Parameters(ParamTemplate):
        TBAS = Float()
        DAGTB = AfgenTrait()
        RVRTB = AfgenTrait()
        BOL50 = Float()
        FALL50 = Float()
        TSOPK = Float() # TSUM until emergence (tsum opkomst)
        TBASE = Float()
        CROP_START_TYPE = Unicode()
        CROP_END_TYPE = Unicode()

    class StateVariables(StatesTemplate):
        DVS = Float()
        BULBSUM = Float()
        BULB = Float()
        EMERGE = Float()
        DOS = Instance(dt.date)
        DOE = Instance(dt.date)
        DOB50 = Instance(dt.date)
        DOF50 = Instance(dt.date)
        STAGE = Enum(["emerging", "vegetative", "reproductive", "mature"])

    class RateVariables(RatesTemplate):
        DTDEV = Float()
        DAYFAC = Float()
        RFR = Float()
        RFRFAC = Float()
        DVR = Float()
        DTSUM = Float()
        DEMERGE = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        if parvalues["CROP_START_TYPE"] == "sowing":
            stage = "emerging"
            dos = day
            doe = None
        else:
            stage = "vegetative"
            dos = None
            doe = day
        self.states = self.StateVariables(kiosk, DVS=0., BULBSUM=0., STAGE=stage,
                                          BULB=0., EMERGE=0., DOS=dos, DOE=doe,
                                          DOB50=None, DOF50=None, publish="DVS", )

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        k = self.kiosk
        p = self.params

        if s.STAGE == "emerging":
            r.DEMERGE = max(0, drv.TEMP - p.TBASE)
            r.DTSUM = 0.
            r.DVR = 0.
        elif s.STAGE == "vegetative":
            r.DTDEV = max(0., drv.TEMP - p.TBAS)
            DL = daylength(day, drv.LAT)
            r.DAYFAC = p.DAGTB(DL)
            r.RFR = exp(-0.222 * k.LAI)
            r.RFRFAC = p.RVRTB(r.RFR)
            r.DEMERGE = 0.
            r.DTSUM = r.DTDEV * r.DAYFAC * r.RFRFAC
            r.DVR = r.DTSUM/p.BOL50
        else:
            r.DEMERGE = 0.
            r.DTDEV = max(0., drv.TEMP - p.TBAS)
            r.DTSUM = r.DTDEV
            r.DVR = r.DTSUM/p.FALL50

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        r = self.rates
        p = self.params

        BULB = 0.
        s.EMERGE += r.DEMERGE * delt
        s.BULBSUM += r.DTSUM * delt
        s.DVS += r.DVR * delt
        if s.STAGE == "emerging":
            if s.EMERGE >= p.TSOPK:
                s.STAGE = "vegetative"
                s.DOE = day
        elif s.STAGE == "vegetative":
            BULB = 0.3 + 100.45 * (exp(-exp(-0.0293*(s.BULBSUM - 91.9))))
            if s.DVS >= 1.0:
                s.STAGE = "reproductive"
                s.DOB50 = day
        elif s.STAGE == "reproductive":
            BULB = 0.3 + 100.45 * (exp(-exp(-0.0293*(s.BULBSUM - 91.9))))
            if s.DVS >= 2.0:
                print("Reached maturity at day %s" % day)
                s.STAGE = "mature"
                s.DOF50 = day
                if p.CROP_END_TYPE == "maturity":
                    self._send_signal(signal=signals.crop_finish, day=day,
                                      finish_type="MATURITY", crop_delete=True)
        else:  # Maturity not more changes in phenological stage
            BULB = 0.3 + 100.45 * (exp(-exp(-0.0293*(s.BULBSUM - 91.9))))

        s.BULB = limit(0., 100., BULB)


class LeafDynamics(SimulationObject):
    """Leaf dynamics for the ALCEPAS crop model.

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

    # parameter for initial LAI as function of plant density
    LAII = Float(-99)
    SLAN = Float(-99)
    SLAR = Float(-99)

    class Parameters(ParamTemplate):
        SLANTB = AfgenTrait()
        SLARTB = AfgenTrait()
        AGECOR = Float(-99)
        METCOR = Float(-99)
        AGEA = Float(-99)
        AGEB = Float(-99)
        AGEC = Float(-99)
        AGED = Float(-99)
        LAGR = Float(-99.)  # Grens tot waar LAI berekend wordt met exp. functie
        GEGR = Float(-99.)  # Totaal droge stof bij LAGR
        LA0 = Float(-99)    # LAI bij opkomst afhankelijk van plantdichtheid (NPL)
        NPL = Float(-99)    # Plant dichtheid
        RGRL = Float(-99)
        GTSLA = Float(-99)
        TTOP = Float(-99)
        TBASE = Float(-99)
        CORFAC = Float(-99)
        TBAS = Float(-99)

    class StateVariables(StatesTemplate):
        LV = Instance(deque)
        SLABC = Instance(deque)
        LVAGE = Instance(deque)
        SPAN = Instance(deque)
        LAIMAX = Float(-99.)
        LAI = Float(-99.)
        WLVG = Float(-99.)
        WLVD = Float(-99.)
        WLV = Float(-99.)
        TSUMEM = Float(-99)

    class RateVariables(RatesTemplate):
        GLV = Float(-99.)
        GLAD = Float(-99.)
        GLA = Float(-99.)
        DLV = Float(-99.)
        SLAT = Float(-99.)
        FYSAGE = Float(-99.)
        SPANT = Float(-99.)
        DTSUMM = Float(-99)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["GLV"])

        # CALCULATE INITIAL STATE VARIABLES
        p = self.params
        self.LAII = p.NPL * p.LA0 * 1.E-4
        self.SLAN = p.SLANTB(p.NPL)
        self.SLAR = p.SLARTB(p.NPL)

        # Initial leaf biomass
        WLVG = 0.
        WLVD = 0.
        WLV = WLVG + WLVD

        # First leaf class (SLA, age and weight)
        LV = deque([WLV])
        LVAGE = deque([0.])
        SLABC = deque([0.])
        SPAN = deque([0.])

        # Initialize StateVariables object
        self.states = self.StateVariables(kiosk, publish=["LAI", "WLV", "WLVG", "WLVD"],
                                          LV=LV, LVAGE=LVAGE, SPAN=SPAN, SLABC=SLABC,
                                          LAIMAX=0., TSUMEM=0., LAI=self.LAII, WLV=WLV,
                                          WLVD=WLVD, WLVG=WLVG)

    def calc_SPAN(self):
        p = self.params
        k = self.kiosk
        SPAN = p.AGECOR * p.METCOR * (p.AGEA + p.AGEB / (1 + p.AGED * k.DVS) + p.AGEC * k.DVS)
        return SPAN

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states
        p = self.params
        k = self.kiosk

        # Growth rate leaves
        # weight of of shoots
        GSH = k.FSH * k.GTW
        # weight of new leaves as fraction of shoots
        r.GLV = k.FLV * GSH

        # Determine how much leaf biomass classes have to die in states.LV,
        # given the a life span > SPAN, these classes will be accumulated
        # in DLV. The dying leaf area is accumulated in GLAD
        # Note that the actual leaf death is imposed on the array LV during the
        # state integration step.
        DLV = 0.0
        GLAD = 0.0
        for lv, lvage, span, sla in zip(s.LV, s.LVAGE, s.SPAN, s.SLABC):
            if lvage > span:
                DLV += lv
                GLAD += sla * lv
        r.DLV = DLV
        r.GLAD = GLAD

        # physiologic ageing of leaves per time step
        DTDEV = max(0, drv.TEMP - p.TBAS)
        r.FYSAGE = DTDEV
        r.SPANT = self.calc_SPAN()

        # Increase and leaf area and SLA
        r.GLA, r.SLAT = self.leaf_area_growth(drv)

    @prepare_states
    def integrate(self, day, delt=1.0):
        params = self.params
        rates = self.rates
        states = self.states

        # --------- leave death ---------
        tLV = array('d', states.LV)
        tSLABC = array('d', states.SLABC)
        tLVAGE = array('d', states.LVAGE)
        tSPAN = array('d', states.SPAN)
        tDLV = rates.DLV

        # leaf death is imposed on leaves by removing leave classes from the
        # right side of the deque.
        for LVweigth in reversed(states.LV):
            if tDLV > 0.:
                if tDLV >= LVweigth:  # remove complete leaf class from deque
                    tDLV -= LVweigth
                    tLV.pop()
                    tLVAGE.pop()
                    tSLABC.pop()
                    tSPAN.pop()
                else:  # Decrease value of oldest (rightmost) leave class
                    tLV[-1] -= tDLV
                    tDLV = 0.
            else:
                break

        # Integration of physiological age
        tLVAGE = deque([age + rates.FYSAGE for age in tLVAGE])
        tLV = deque(tLV)
        tSLABC = deque(tSLABC)
        tSPAN = deque(tSPAN)

        # --------- leave growth ---------
        # new leaves in class 1
        tLV.appendleft(rates.GLV)
        tSLABC.appendleft(rates.SLAT)
        tLVAGE.appendleft(0.)
        tSPAN.appendleft(rates.SPANT)

        # calculation of new leaf area
#        states.LAI = sum([lv * sla for lv, sla in zip(tLV, tSLABC)])
        states.LAI += rates.GLA
        states.LAIMAX = max(states.LAI, states.LAIMAX)

        # Update leaf biomass states
        states.WLVG = sum(tLV)
        states.WLVD += rates.DLV
        states.WLV = states.WLVG + states.WLVD

        # Store final leaf biomass deques
        self.states.LV = tLV
        self.states.SLABC = tSLABC
        self.states.LVAGE = tLVAGE
        self.states.SPAN = tSPAN

        self.states.TSUMEM += self.rates.DTSUMM

    def leaf_area_growth(self, drv):
        # Computes daily increase of leaf area index (ha leaf/ ha ground/ d)
        p = self.params
        k = self.kiosk
        s = self.states
        r = self.rates

        DTEFF = limit(p.TBASE, p.TTOP, drv.TEMP)
        DTR = drv.IRRAD
        if DTEFF > 0:
            r.DTSUMM = 1./((1./DTEFF) + p.CORFAC*(1/(0.5*0.000001*DTR)))
        else:
            r.DTSUMM = 0.

        if s.LAI < p.LAGR and k.TADRW < p.GEGR:
            # leaf growth during juvenile growth:
            SLA = (self.SLAN + self.SLAR * p.GTSLA ** k.DVS) * 1/100000.
            GLA = self.LAII * p.RGRL * r.DTSUMM * exp(p.RGRL * s.TSUMEM)
            # Adjust SLA for youngest leaves under exponential growth conditions
            if r.GLV > 0.:
                SLA = GLA/r.GLV
        else:
            # leaf growth during mature plant growth:
            SLA = (self.SLAN + self.SLAR * p.GTSLA ** k.DVS) * 1/100000.
            GLA = (SLA * r.GLV)

        # correct for leaf death
        GLA = GLA - r.GLAD

        return GLA, SLA


class RootDynamics(SimulationObject):

    class RateVariables(RatesTemplate):
        GRT = Float()

    class StateVariables(StatesTemplate):
        WRT = Float()

    def initialize(self, day, kiosk, parvalues):
        self.rates = self.RateVariables(kiosk, publish=["GRT"])
        self.states = self.StateVariables(kiosk, WRT=0., publish=["WRT"])

    @prepare_rates
    def calc_rates(self, day, drv):
        k = self.kiosk
        r = self.rates
        r.GRT = k.FRT * k.GTW

    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        s = self.states
        s.WRT += r.GRT * delt


class BulbDynamics(SimulationObject):

    class RateVariables(RatesTemplate):
        GSO = Float()

    class StateVariables(StatesTemplate):
        WSO = Float()

    def initialize(self, day, kiosk, parvalues):
        self.rates = self.RateVariables(kiosk, publish=["GSO"])
        self.states = self.StateVariables(kiosk, WSO=0., publish=["WSO"])

    @prepare_rates
    def calc_rates(self, day, drv):
        k = self.kiosk
        r = self.rates
        # increase in weight of shoots
        GSH = k.FSH * k.GTW
        # Increase in weight of bulb
        r.GSO = k.FSO * GSH

    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        s = self.states
        s.WSO += r.GSO * delt


class ALCEPAS(SimulationObject):
    leafdynamics = Instance(SimulationObject)
    phenology = Instance(SimulationObject)
    partitioning = Instance(SimulationObject)
    assimilation = Instance(SimulationObject)
    respiration = Instance(SimulationObject)
    rootdynamics = Instance(SimulationObject)
    bulbdynamics = Instance(SimulationObject)

    class Parameters(ParamTemplate):
        ASRQSO = Float()
        ASRQRT = Float()
        ASRQLV = Float()

    class RateVariables(RatesTemplate):
        GTW = Float()

    class StateVariables(StatesTemplate):
        TADRW = Float()

    def initialize(self, day, kiosk, parvalues):
        self.leafdynamics = LeafDynamics(day, kiosk, parvalues)
        self.phenology = Phenology(day, kiosk, parvalues)
        self.partitioning = Partitioning(day, kiosk, parvalues)
        self.assimilation = Assimilation(day, kiosk, parvalues)
        self.respiration = Respiration(day, kiosk, parvalues)
        self.rootdynamics = RootDynamics(day, kiosk, parvalues)
        self.bulbdynamics = BulbDynamics(day, kiosk, parvalues)

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["GTW"])
        self.states = self.StateVariables(kiosk, TADRW=0., publish=["TADRW"])

    @prepare_rates
    def calc_rates(self, day, drv):
        p = self.params
        k = self.kiosk
        r = self.rates

        # phenological development
        self.phenology.calc_rates(day, drv)
        if self.get_variable("STAGE") == "emerging":
            self.touch()  # No need to continue before emergence
            return

        # Assimilation and respiration
        GPHOT = self.assimilation(day, drv)
        MAINT = self.respiration(day, drv)

        # Partitioning and assimilate requirements for dry matter conversion (kgCH20 / kgDM)
        self.partitioning(day, drv)
        ASRQ = k.FSH * (p.ASRQLV * k.FLV + p.ASRQSO * k.FSO) + p.ASRQRT * k.FRT
        # Total dry matter growth
        r.GTW = (GPHOT - MAINT) / ASRQ

        # Partitioning and dynamics of different plant organs
        self.leafdynamics.calc_rates(day, drv)
        self.rootdynamics.calc_rates(day, drv)
        self.bulbdynamics.calc_rates(day, drv)

        self.check_carbon_balance(day)

    def check_carbon_balance(self, day):
        r = self.rates
        k = self.kiosk
        if r.GTW - (k.GSO + k.GLV + k.GRT) > 0.0001:
            raise exc.CarbonBalanceError("Carbon balance not closing on day %s" % day)

    @prepare_states
    def integrate(self, day, delt=1.0):
        k = self.kiosk
        s = self.states

        self.leafdynamics.integrate(day, delt)
        self.phenology.integrate(day, delt)
        self.rootdynamics.integrate(day, delt)
        self.bulbdynamics.integrate(day, delt)

        s.TADRW = k.WLV + k.WSO + k.WRT
