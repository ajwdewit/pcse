# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
"""
Class to calculate various nutrient relates stress factors:
    NNI      nitrogen nutrition index   
    NPKI     NPK nutrition index (= minimum of N/P/K-index)
    NPKREF   assimilation reduction factor based on NPKI
"""

from ...traitlets import Float
from ...util import limit, AfgenTrait
from ...base import ParamTemplate, SimulationObject, RatesTemplate
from ...decorators import prepare_rates

class NPK_Stress(SimulationObject):
    """Implementation of N stress calculation through [N]nutrition index.

    Stress factors are calculated based on the mass concentration of N in
    the leaf and stem biomass of the plant. For each pool of nutrients, four
    concentrations are calculated based on the biomass for leaves and stems:
    - the actual concentration based on the actual amount of nutrients
      divided by the actual leaf and stem biomass.
    - The maximum concentration, being the maximum that the plant can absorb
      into its leaves and stems.
    - The critical concentration, being the concentration that is needed to
      maintain growth rates that are not limited by N.For N, the critical
      concentration can be lower than the maximum concentration. This 
      concentration is sometimes called 'optimal concentration'.
    - The residual concentration which is the amount that is locked
      into the plant structural biomass and cannot be mobilized anymore.

    The stress index (SI) is determined as a simple ratio between those
    concentrations according to:

    :math:`SI = (C_{a} - C_{r})/(C_{c} - C_{r})`

    with subscript `a`, `r` and `c` being the actual, residual and critical
    concentration for the N.
    This equation results in the nitrogen nutrition index (NNI). Finally, the reduction 
    factor for assimilation (NPKREF) is calculated using the reduction factor for
    light use efficiency (NLUE_NPK).

    **Simulation parameters**

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
    NNI       Nitrogen nutrition index                          Y     -
    RFNPK     Reduction factor for |CO2| assimlation            N     -
              based on NPKI and the parameter NLUE_NPK
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

        NCRIT_FR = Float(-99.)   # optimal N concentration as fraction of maximum N concentration
        NMAXRT_FR = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        NMAXST_FR = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NLUE_NPK = Float(-99.)  # coefficient for the reduction of RUE due to nutrient (N-P-K) stress
        NMAXSO = Float(-99.)
        RGRLAI_MIN = Float(-99.)
        RGRLAI = Float(-99.)

    class RateVariables(RatesTemplate):
        NNI = Float()
        NPKI = Float()
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
        self.rates = self.RateVariables(kiosk, publish=["NPKI", "NNI", "NSLLV", "RFRGRL"])

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
        
        # Total vegetative living above-ground biomass (kg DM ha-1)
        VBM = k.WLV + k.WST
      
        # Critical (Optimal) NPK amount in vegetative above-ground living biomass
        # and its NPK concentration
        NcriticalLV  = p.NCRIT_FR * NMAXLV * k.WLV
        NcriticalST  = p.NCRIT_FR * NMAXST * k.WST
        
        # if above-ground living biomass = 0 then optimum = 0
        if VBM > 0.:
            NcriticalVBM = (NcriticalLV + NcriticalST)/VBM
        else:
            NcriticalVBM = 0.

        # NPK concentration in total vegetative living per kg above-ground
        # biomass  (kg N/P/K kg-1 DM)
        # if above-ground living biomass = 0 then concentration = 0
        if VBM > 0.:
            NconcentrationVBM  = (k.NamountLV + k.NamountST)/VBM
        else:
            NconcentrationVBM = 0.

        # Residual NPK concentration in total vegetative living above-ground
        # biomass  (kg N kg-1 DM)
        # if above-ground living biomass = 0 then residual concentration = 0
        if VBM > 0.:
            NresidualVBM = (k.WLV * p.NRESIDLV + k.WST * p.NRESIDST)/VBM
        else:
            NresidualVBM = PresidualVBM = KresidualVBM = 0.
            
        if (NcriticalVBM - NresidualVBM) > 0.:
            r.NNI = limit(0.001, 1.0, (NconcentrationVBM - NresidualVBM)/(NcriticalVBM - NresidualVBM))
        else:
            r.NNI = 0.001
            
        r.NPKI = r.NNI

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

        # Nutrient reduction factor for assimilation
        r.RFNPK = limit(0., 1.0, 1. - (p.NLUE_NPK * (1.0001 - r.NPKI) ** 2))
         
        return r.NNI, r.NPKI, r.RFNPK
