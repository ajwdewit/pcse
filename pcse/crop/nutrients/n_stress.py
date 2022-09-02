# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
"""
Class to calculate nitrogen stress factors:

"""

from ...traitlets import Float
from ...util import limit, AfgenTrait
from ...base import ParamTemplate, SimulationObject, RatesTemplate
from ...decorators import prepare_rates

class N_Stress(SimulationObject):
    """Implementation of N stress calculation through [N]nutrition index.

    HB 20220405 A lot of changes have been done in this subroutine. It needs to be redocumented.

    ============  ============================================= ======================
     Name          Description                                   Unit
    ============  ============================================= ======================
    NMAXLV_TB      Maximum N concentration in leaves as         kg N kg-1 dry biomass
                   function of DVS
    NMAXRT_FR      Maximum N concentration in roots as fraction -
                   of maximum N concentration in leaves
    NMAXST_FR      Maximum N concentration in stems as fraction -
                   of maximum N concentration in leaves
    NCRIT_FR       Critical N concentration as fraction of      -
                   maximum N concentration for vegetative
                   plant organs as a whole (leaves + stems)
    NRESIDLV       Residual N fraction in leaves                kg N kg-1 dry biomass
    NRESIDST       Residual N fraction in stems                 kg N kg-1 dry biomass

    NLUE_NPK       Coefficient for the reduction of RUE due     -
                   to nutrient (N-P-K) stress
    ============  ============================================= ======================

    **Rate variables**

    The rate variables here are not real rate variables in the sense that they are derived
    state variables and do not represent a rate. However, as they are directly used
    in the rate variable calculation it is logical to put them here.

    =======  ================================================= ==== ==============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ==============
    NSTRESS   Ratio of maximum to actual aboveground crop N     Y     -
    =======  ================================================= ==== ==============


    **External dependencies:**

    ==========  =================================== =====================  ==============
     Name        Description                         Provided by            Unit
    ==========  =================================== =====================  ==============
    DVS          Crop development stage              DVS_Phenology           -
    WST          Dry weight of living stems          WOFOST_Stem_Dynamics  |kg ha-1|
    WLV          Dry weight of living leaves         WOFOST_Leaf_Dynamics  |kg ha-1|
    NamountLV    Amount of N in leaves               NPK_Crop_Dynamics     |kg ha-1|
    NamountST    Amount of N in stems                NPK_Crop_Dynamics     |kg ha-1|
    ==========  =================================== =====================  ==============
    """

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs
        NSLLV_TB = AfgenTrait()      # N stress multiplication factor for leaf death
        NMAXRT_FR = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        NMAXST_FR = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NMAXSO = Float(-99.)
        RGRLAI_MIN = Float(-99.)
        RGRLAI = Float(-99.)

    class RateVariables(RatesTemplate):
        RFNPK = Float()
        NSLLV = Float()
        RFRGRL = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: current date
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: ParameterProvider with parameter key/value pairs
        """

        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish = ["NSLLV", "RFRGRL"])

    @prepare_rates
    def __call__(self, day, drv):
        """

        :param day: the current date
        :param drv: the driving variables
        :return: A tuple (NNI, NPKI, NPKREF)
        """
        p = self.params
        r = self.rates
        k = self.kiosk

        # Maximum NPK concentrations in leaves (kg N kg-1 DM)
        NMAXLV = p.NMAXLV_TB(k.DVS)

        # Maximum NPK concentrations in stems (kg N kg-1 DM)
        NMAXST = p.NMAXST_FR * NMAXLV
        
        # Calculate multiplication factor of leaf death due to N stress
        NamountABG = k.NamountLV + k.NamountST + k.NamountSO
        NamountABGMX = k.WLV * NMAXLV + k.WST * NMAXST + k.WSO * p.NMAXSO

        if NamountABGMX / NamountABG <= 1:
            NstressIndexDLV = 1.
        elif NamountABGMX / NamountABG > 2:
            NstressIndexDLV = 2.
        else:
            NstressIndexDLV = NamountABGMX / NamountABG 
        
        r.NSLLV = p.NSLLV_TB(NstressIndexDLV)

        # Calculate reduction factor of leaf growth rate in exponential growth phase
        if(k.WLV > 0):
            NconcentrationLV = k.NamountLV / k.WLV
        else:
            NconcentrationLV = 0.

        NstressIndexRGRLAI = max(0, min(1, (NconcentrationLV - 0.9 * NMAXLV) / (NMAXLV - 0.9 * NMAXLV)))
        r.RFRGRL = 1 - (1.-NstressIndexRGRLAI)*(p.RGRLAI-p.RGRLAI_MIN) / p.RGRLAI