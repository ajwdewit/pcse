# -*- coding: utf-8 -*-
# Copyright (c) 2021 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2021
"""Implementation of the LINGRA grassland simulation model

This module provides an implementation of the LINGRA (LINtul GRAssland)
simulation model for grasslands as described by Schapendonk et al. 1998
(https://doi.org/10.1016/S1161-0301(98)00027-6) for use within the
Python Crop Simulation Environment.
"""
from math import exp, log

from pcse.base import SimulationObject, ParamTemplate, StatesTemplate, RatesTemplate
from pcse.traitlets import Float, List, Bool, Instance, Integer
from pcse.util import AfgenTrait, limit
from pcse.decorators import prepare_states, prepare_rates
from pcse.crop.evapotranspiration import Evapotranspiration
from pcse.crop.root_dynamics import Simple_Root_Dynamics
import pcse.signals


class SourceLimitedGrowth(SimulationObject):
    """Calculates the source-limited growth rate for grassland based on radiation and
    temperature as driving variables and possibly limited by soil moisture or
    leaf nitrogen content.The latter is based on static values for current and
    maximum N concentrations and is mainly there for connecting an N module in the
    future.

    This routine uses a light use efficiency (LUE) approach where the LUE is adjusted
    for effects of temperature and radiation level. The former is need as photosynthesis
    has a clear temperature response. The latter is required as photosynthesis rate
    flattens off at higher radiation levels which leads to a lower 'apparent' light use \
    efficiency. The parameter `LUEreductionRadiationTB` is a crude empirical correction
    for this effect.

    Note that a reduction in growth rate due to soil moisture is obtained through the
    reduction factor for transpiration (RFTRA).

    This module does not provide any true rate variables, but returns the computed
    growth rate directly to the calling routine through __call__().

    *Simulation parameters*:

    =======================  =============================================  ==============
     Name                      Description                                     Unit
    =======================  =============================================  ==============
    KDIFTB                    Extinction coefficient for diffuse visible        -
                              as function of DVS.
    CO2A                      Atmospheric CO2 concentration                    ppm
    LUEreductionSoilTempTB    Reduction function for light use efficiency      C, -
                              as a function of soil temperature.
    LUEreductionRadiationTB   Reduction function for light use efficiency      MJ, -
                              as a function of radiation level.
    LUEmax                    Maximum light use efficiency.
    =======================  =============================================  ==============


    *Rate variables*

    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
    RF_RadiationLevel     Reduction factor for light use efficiency       -
                          due to the radiation level
    RF_RadiationLevel     Reduction factor for light use efficiency       -
                          due to the radiation level
    LUEact                The actual light use efficiency                g /(MJ PAR)
    ===================  =============================================  ===============

    *Signals send or handled*

    None

    *External dependencies:*

    ===============  =================================== ==============================
     Name             Description                         Provided by
    ===============  =================================== ==============================
    DVS               Crop development stage              pylingra.LINGRA
    TemperatureSoil   Soil Temperature                    pylingra.SoilTemperature
    RFTRA             Reduction factor for transpiration  pcse.crop.Evapotranspiration
    ===============  =================================== ==============================
    """

    class Parameters(ParamTemplate):
        KDIFTB = AfgenTrait()
        LUEreductionSoilTempTB = AfgenTrait()
        LUEreductionRadiationTB = AfgenTrait()
        CO2A = Float()
        LUEmax = Float()

    class RateVariables(RatesTemplate):
        RF_Temperature = Float()
        RF_RadiationLevel = Float()
        LUEact = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(self.kiosk, publish=["RF_Temperature"])

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk

        # From J/m2/d to MJ/m2/d
        DTR = drv.IRRAD / 1.E+6
        PAR = DTR * 0.50

        # Photosynthesis reduction factors for temperature and radiation level
        r.RF_Temperature = p.LUEreductionSoilTempTB(k.TemperatureSoil)
        r.RF_RadiationLevel = p.LUEreductionRadiationTB(DTR)

        # Fraction of light interception
        FINT = (1.-exp(-p.KDIFTB(k.DVS) * k.LAI))

        # Total intercepted photosynthetically active
        # radiation, MJ m-2 d-1
        PARINT = FINT * PAR

        # Light use efficiency corrected for temp and radiation level, g MJ PAR-1
        LUEpot = p.LUEmax * r.RF_Temperature * r.RF_RadiationLevel

        # LUE corrected for transpiration stress
        r.LUEact = LUEpot * k.RFTRA

        if k.dWeightHARV == 0.:  # No grass harvest today, normal growth
            # (10: conversion from g m-2 d-1 to kg ha-1 d-1)
            GrowthSource = r.LUEact * PARINT * (1. + 0.8 * log(p.CO2A / 360.)) * 10.
        else:
            # growth is zero at day when mowing occurs
            GrowthSource = 0.

        return GrowthSource


class SinkLimitedGrowth(SimulationObject):
    """Calculates the sink-limited growth rate for grassland assuming a temperature
    driven maximum leaf elongation rate multiplied by the number of tillers. The
    conversion to growth in kg/ha dry matter is done by dividing by the specific
    leaf area (SLA).

    Besides the sink-limited growth rate, this class also computes the change
    in tiller number taking into account the growth rate, death rate and number
    of days after defoliation due to harvest.

    *Simulation parameters*:

    =======================  =============================================  ==============
     Name                      Description                                     Unit
    =======================  =============================================  ==============
    TempBase                  Base temperature for leaf development and
                              grass phenology                                  C
    LAICrit                   Cricical leaf area beyond which leaf death
                              due to self-shading occurs                       -
    SiteFillingMax            Maximum site filling for new buds             tiller/leaf-1
    SLA                       Specific leaf area                             ha/kg
    TSUMmax                   Temperature sum to max development stage        C.d
    TillerFormRateA0          A parameter in the equation for tiller
                              formation rate valid up till 7 days after
                              harvest
    TillerFormRateB0          B parameter in the equation for tiller
                              formation rate valid up till 7 days after
                              harvest
    TillerFormRateA8          A parameter in the equation for tiller
                              formation rate starting from 8 days after
                              harvest
    TillerFormRateB8          B parameter in the equation for tiller
                              formation rate starting from 8 days after
                              harvest
    =======================  =============================================  ==============


    *Rate variables*

    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
    dTillerNumber         Change in tiller number                        tillers/m2/d
                          due to the radiation level
    dLeafLengthPot        Potential change in leaf length. Later on      cm/d
                          the actual change in leaf length will be
                          computed taking source limitation into
                          account.
    LAIGrowthSink         Growth of LAI based on sink-limited growth     d-1
                          rate.
    ===================  =============================================  ===============

    *Signals send or handled*

    None

    *External dependencies:*

    ===============  =================================== ==============================
     Name             Description                         Provided by
    ===============  =================================== ==============================
    DVS               Crop development stage              pylingra.LINGRA
    LAI               Leaf Area Index                     pylingra.LINGRA
    TemperatureSoil   Soil Temperature                    pylingra.SoilTemperature
    RF_Temperature    Reduction factor for LUE based on
                      temperature                         pylingra.SourceLimitedGrowth
    TillerNumber      Actual number of tillers            pylingra.LINGRA
    LVfraction        Fraction of assimilates going to    pylingra.LINGRA
                      leaves
    dWeightHARV       Change in harvested weight          pylingra.LINGRA
                      (indicates that a harvest took
                      place today)
    ===============  =================================== ==============================
"""

    class Parameters(ParamTemplate):
        TempBase = Float()
        LAIcrit = Float()
        SiteFillingMax = Float()
        SLA = Float()
        TillerFormRateA0 = Float()
        TillerFormRateB0 = Float()
        TillerFormRateA8 = Float()
        TillerFormRateB8 = Float()
        TSUMmax = Float()

    class RateVariables(RatesTemplate):
        dTillerNumber = Float()
        dLeafLengthPot = Float()
        LAIGrowthSink = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["dTillerNumber", "dLeafLengthPot"])

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk

        # Temperature dependent leaf appearance rate, according to
        # (Davies and Thomas, 1983), soil temperature (TemperatureSoil) is used as
        # driving force which is estimated from a 10 day running average
        LeafAppRate = k.TemperatureSoil * 0.01 if k.RF_Temperature > 0. else 0.
        # derive tiller rate
        r.dTillerNumber = self._calc_tillering_rate(LeafAppRate)

        # Leaf elongation rate affected by temperature: cm day-1 tiller-1
        r.dLeafLengthPot = 0.83 * log(max(drv.TEMP, 2.)) - 0.8924 if (drv.TEMP - p.TempBase) > 0. else 0.

        # Rate of sink limited leaf growth, unit of TillerNumber is tillers m-2
        # 1.0E-8 is conversion from cm-2 to ha-1, ha leaf ha ground-1 d-1
        r.LAIgrowthSink = (k.TillerNumber * 1.0E4 * (r.dLeafLengthPot * 0.3)) * 1.0E-8
        # Conversion of leaf growth rate to total sink limited carbon demand using SLA
        # in kg leaf ha-1 d-1
        GrowthSink = r.LAIgrowthSink * (1./p.SLA) * (1./k.LVfraction) if k.dWeightHARV <= 0. else 0.

        return GrowthSink

    def _calc_tillering_rate(self, LeafAppRate):
        k = self.kiosk
        p = self.params
        # Actual site filling equals maximum site filling without N stress
        SiteFillingAct =  p.SiteFillingMax

        if k.DaysAfterHarvest < 8.:
            # Relative rate of tiller formation when defoliation less
            # than 8 days ago, tiller tiller-1 d-1
            TillerFormationRate = max(0., p.TillerFormRateA0 - p.TillerFormRateB0 * k.LAI) * k.RF_Temperature
        else:
            # Relative rate of tiller formation when defoliation is more
            # than 8 days ago, tiller tiller-1 d-1
            TillerFormationRate = limit(0., SiteFillingAct, p.TillerFormRateA8 - p.TillerFormRateB8 * k.LAI) * k.RF_Temperature

        # Relative death rate of tillers due to self-shading (DTILD), tiller tiller-1 d-1
        TillerDeathRate = max(0.01 * (1. + k.TSUM / p.TSUMmax), 0.05 * (k.LAI - p.LAIcrit) / p.LAIcrit)

        # Change in Tiller number
        if k.TillerNumber <= 14000.:
            dTillerNumber = (TillerFormationRate - TillerDeathRate) * LeafAppRate * k.TillerNumber
        else:
            dTillerNumber = -TillerDeathRate * LeafAppRate * k.TillerNumber

        return dTillerNumber


class SoilTemperature(SimulationObject):
    """Calculates the soil temperature in the upper 10 cm using a 10-day
    moving average of the daily average air temperature at 2m.

    *Simulation parameters*:

    =======================  =============================================  ==============
     Name                      Description                                     Unit
    =======================  =============================================  ==============
    SoilTemperatureInit       Initial soil temperature                         C
    =======================  =============================================  ==============


    *Rate variables*
    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
    dTemperatureSoil      Change in soil temperature                      C/d
    ===================  =============================================  ===============

    *State variables*
    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
    TemperatureSoil       Actual soil temperature                         C
    ===================  =============================================  ===============

    *Signals send or handled*

    None

    *External dependencies:*

    None
    """

    class Parameters(ParamTemplate):
        TemperatureSoilinit = Float()

    class StateVariables(StatesTemplate):
        TemperatureSoil = Float()

    class RateVariables(RatesTemplate):
        dTemperatureSoil = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk)
        self.states = self.StateVariables(kiosk, TemperatureSoil=self.params.TemperatureSoilinit,
                                          publish=["TemperatureSoil"])

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states

        # soil temperature changes
        r.dTemperatureSoil = (drv.TEMP - s.TemperatureSoil) / 10.

    @prepare_states
    def integrate(self, day, delt=1.0):
        self.states.TemperatureSoil += self.rates.dTemperatureSoil * delt


class LINGRA(SimulationObject):
    """Top level implementation of LINGRA, integrating all components

    This class integrates all components from the LINGRA model and includes the
    main state variables related to weights of the different biomass pools, the
    leaf area, tiller number and leaf length. The integrated components include the
    implementations for source/sink limited growth, soil temperature,
    evapotranspiration and root dynamics. The latter two are taken from WOFOST in
    order to avoid duplication of code.

    Compared to the original code from Schapendonk et al. (1998) several improvements
    have been made:

    - an overall restructuring of the code, removing unneeded variables and renaming
      the remaining variables to have more readable names.
    - A clearer implementation of sink/source limited growth including the use of
      reserves
    - the potential leaf elongation rate as calculated by the Sink-limited growth
      module is now corrected for actual growth. Thereby avoiding unlimited leaf
      growth under water-stressed conditions which led to unrealistic results.

    *Simulation parameters*:

    =======================  =============================================  ==============
     Name                      Description                                     Unit
    =======================  =============================================  ==============
     LAIinit                  Initial leaf area index                           -
     TillerNumberinit         Initial number of tillers                      tillers/m2
     WeightREinit             Initial weight of reserves                     kg/ha
     WeightRTinit             Initial weight of roots                        kg/ha
     LAIcrit                  Critical LAI for death due to self-shading     -
     RDRbase                  Background relative death rate for roots       d-1
     RDRShading               Max relative death rate of leaves due to       d-1
                              self-shading
     RDRdrought               Max relative death rate of leaves due to
                              drought stress                                 d-1
     SLA                      Specific leaf area                             ha/kg
     TempBase                 Base temperature for photosynthesis and        C
                              development
     PartitioningRootsTB      Partitioning fraction to roots as a            -, -
                              function of the reduction factor for
                              transpiration (RFTRA)
     TSUMmax                  Temperature sum to max development stage       C.d
    =======================  =============================================  ==============


    *Rate variables*

    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
    dTSUM                 Change in temperature sum for development       C
    dLAI                  Net change in Leaf Area Index                   d-1
    dDaysAfterHarvest     Change in Days after Harvest                    -
    dCuttingNumber        Change in number of cuttings (harvests)         -
    dWeightLV             Net change in leaf weight                       kg/ha/d
    dWeightRE             Net change in reserve pool                      kg/ha/d
    dLeafLengthAct        Change in actual leaf length                    cm/d
    LVdeath               Leaf death rate                                 kg/ha/d
    LVgrowth              Leaf growth rate                                kg/ha/d
    dWeightHARV           Change in harvested dry matter                  kg/ha/d
    dWeightRT             Net change in root weight                       kg/ha/d
    LVfraction            Fraction partitioned to leaves                  -
    RTfraction            Fraction partitioned to roots                   -
    ===================  =============================================  ===============

    *State variables*

    ===================  =============================================  ===============
     Name                 Description                                     Unit
    ===================  =============================================  ===============
     TSUM                 Temperature sum                                  C d
     LAI                  Leaf area Index                                  -
     DaysAfterHarvest     number of days after harvest                     d
     CuttingNumber        number of cuttings (harvests)                    -
     TillerNumber         Tiller number                                    tillers/m2
     WeightLVgreen        Weight of green leaves                           kg/ha
     WeightLVdead         Weight of dead leaves                            kg/ha
     WeightHARV           Weight of harvested dry matter                   kg/ha
     WeightRE             Weight of reserves                               kg/ha
     WeightRT             Weight of roots                                  kg/ha
     LeafLength           Length of leaves                                 kg/ha
     WeightABG            Total aboveground weight (harvested +            kg/ha
                          current)
     SLAINT               Integrated SLA during the season                 ha/kg
     DVS                  Development stage                                -
    ===================  =============================================  ===============

    *Signals sent or handled*

    Mowing of grass will take place when a `pcse.signals.mowing` event is broadcasted.
    This will reduce the amount of living leaf weight assuming that a certain
    amount of biomass will remain on the field (this is a parameter on the MOWING
    event).

    *External dependencies:*

    ===============  =================================== ====================================
     Name             Description                         Provided by
    ===============  =================================== ====================================
    RFTRA             Reduction factor for transpiration  pcse.crop.Evapotranspiration
    dLeafLengthPot    Potential growth in leaf length     pcse.crop.lingra.SinkLimitedGrowth
    dTillerNumber     Change in tiller number             pcse.crop.lingra.SinkLimitedGrowth
    ===============  =================================== ====================================
    """

    WeightLV_remaining = Float()
    _flag_MOWING = Bool(False)

    source_limited_growth = Instance(SimulationObject)
    sink_limited_growth = Instance(SimulationObject)
    soil_temperature = Instance(SimulationObject)
    evapotranspiration = Instance(SimulationObject)
    root_dynamics = Instance(SimulationObject)

    class Parameters(ParamTemplate):
        LAIinit = Float()
        TillerNumberinit = Float()
        WeightREinit = Float()
        WeightRTinit = Float()
        LAIcrit = Float()
        RDRbase = Float()
        SLA = Float()
        TempBase = Float()
        PartitioningRootsTB = AfgenTrait()
        TSUMmax = Float()
        RDRshading = Float()
        RDRdrought = Float()

    class RateVariables(RatesTemplate):
        dTSUM = Float()
        dLAI = Float()
        dDaysAfterHarvest = Integer()
        dCuttingNumber = Integer()
        dWeightLV = Float()
        dWeightRE = Float()
        dLeafLengthAct = Float()
        LVdeath = Float()
        LVgrowth = Float()
        dWeightHARV = Float()
        dWeightRT = Float()
        LVfraction = Float()
        RTfraction = Float()

    class StateVariables(StatesTemplate):
        TSUM = Float()
        LAI = Float()
        DaysAfterHarvest = Integer()
        CuttingNumber = Integer()
        TillerNumber = Float()
        WeightLVgreen = Float()
        WeightLVdead = Float()
        WeightHARV = Float()
        WeightRE = Float()
        WeightRT = Float()
        LeafLength = Float()
        WeightABG = Float()
        SLAINT = Float()
        DVS = Float()

    def initialize(self, day, kiosk, parvalues):

        self.source_limited_growth = SourceLimitedGrowth(day, kiosk, parvalues)
        self.sink_limited_growth = SinkLimitedGrowth(day, kiosk, parvalues)
        self.soil_temperature = SoilTemperature(day, kiosk, parvalues)
        self.evapotranspiration = Evapotranspiration(day, kiosk, parvalues)
        self.root_dynamics = Simple_Root_Dynamics(day, kiosk, parvalues)

        self.params = self.Parameters(parvalues)
        p = self.params
        s = {"TSUM": 0.,
             "LAI": p.LAIinit,
             "DaysAfterHarvest": 0,
             "CuttingNumber": 0,
             "TillerNumber": p.TillerNumberinit,
             "WeightLVgreen": p.LAIinit / p.SLA,
             "WeightLVdead": 0.,
             "WeightHARV": 0.,
             "WeightRE": p.WeightREinit,
             "WeightRT": p.WeightRTinit,
             "LeafLength": 0.,
             "WeightABG": p.LAIinit / p.SLA,
             "SLAINT": p.SLA,
             "DVS": 0.0}
        pub = ["LAI", "WeightRE", "LeafLength", "TillerNumber", "TSUM",
               "DaysAfterHarvest", "DVS"]
        self.states = self.StateVariables(kiosk, **s, publish=pub)
        self.rates = self.RateVariables(kiosk, publish=["dWeightHARV", "LVfraction", "RTfraction"])
        self._connect_signal(self._on_MOWING, signal=pcse.signals.mowing)

    @prepare_rates
    def calc_rates(self, day, drv):
        p = self.params
        r = self.rates
        k = self.kiosk
        s = self.states

        r.dTSUM = max(drv.TEMP - p.TempBase, 0.)

        self.soil_temperature.calc_rates(day, drv)
        self.root_dynamics.calc_rates(day, drv)
        self.evapotranspiration(day, drv)

        # grassland management options for mowing
        if self._flag_MOWING:
            r.dWeightHARV = max(0, s.WeightLVgreen - self.WeightLV_remaining)
            r.dDaysAfterHarvest = -s.DaysAfterHarvest
            r.dCuttingNumber = 1
        else:
            r.dCuttingNumber = 0
            r.dWeightHARV = 0.
            if s.CuttingNumber > 0 or r.dCuttingNumber == 1:
                r.dDaysAfterHarvest = 1
            else:
                r.dDaysAfterHarvest = 0

        # *** Death rates leaves ***
        # Relative death rate of leaves due to drought stress, d-1
        RDRdrought = limit(0., p.RDRdrought, p.RDRdrought * (1.-k.RFTRA))
        # Relative death rate of leaves due to self-shading, d-1
        RDRshading = limit(0., p.RDRshading, p.RDRshading * (s.LAI-p.LAIcrit)/p.LAIcrit)
        # Actual relative death rate of leaves is sum of base death
        # rate plus maximum of death rates of shading and drought, d-1
        RDRtotal = p.RDRbase + max(RDRshading, RDRdrought)
        # Actual death rate of leaf area, due to relative death
        # rate of leaf area or rate of change due to cutting, ha ha-1 d-1
        if self._flag_MOWING:
            LAIdeath = r.dWeightHARV * s.SLAINT
        else:
            LAIdeath = s.LAI * (1. - exp(-RDRtotal))

        # Fraction of dry matter allocated to roots/leaves, kg kg-1
        r.RTfraction = p.PartitioningRootsTB(k.RFTRA)
        r.LVfraction = 1. - r.RTfraction

        # *** Growth rates ***
        GrowthSource = self.source_limited_growth(day, drv)
        GrowthSink = self.sink_limited_growth(day, drv)

        # Actual growth switches between sink- and source limitation
        if GrowthSource < GrowthSink:  # source limited growth
            gap = GrowthSink - GrowthSource  # gap in assimilates
            dWeightRE = min(s.WeightRE, gap)  # available reserves
            GrowthAct = GrowthSource + dWeightRE
            r.dWeightRE = -dWeightRE
        else:  # Sink limited growth
            r.dWeightRE = GrowthSource - GrowthSink  # surplus in assimilates
            GrowthAct = GrowthSink

        # Actual growth rate of leaf area, ha ha-1 d-1
        LAIgrowthAct = GrowthAct * r.LVfraction * p.SLA

        # rate of change of dry weight of green leaves due to
        # growth and senescence of leaves or periodical harvest, kg ha-1
        r.LVgrowth = GrowthAct * r.LVfraction if r.dWeightHARV <= 0 else 0.

        # Actual death rate of leaf biomass, kg ha-1 d-1 incl. harvested leaves
        r.LVdeath = LAIdeath / s.SLAINT

        # Change in LAI
        r.dLAI = LAIgrowthAct - LAIdeath

        # Change in green leaf weight
        r.dWeightLV = r.LVgrowth - r.LVdeath

        # Actual growth rate of roots, kg ha-1 d-1
        r.dWeightRT = GrowthAct * r.RTfraction

        # Actual change in leaf length
        if r.dWeightHARV > 0:  # Correction for harvesting
            r.dLeafLengthAct = -s.LeafLength
        else:
            if GrowthSink > 0:
                # Correction for source limitation
                r.dLeafLengthAct = k.dLeafLengthPot * GrowthAct/GrowthSink
            else:
                r.dLeafLengthAct = 0

        self._flag_MOWING = False

    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        k = self.kiosk
        s = self.states
        p = self.params

        s.TSUM += r.dTSUM * delt
        s.LAI += r.dLAI * delt
        s.DaysAfterHarvest += r.dDaysAfterHarvest
        s.CuttingNumber += r.dCuttingNumber
        s.TillerNumber += k.dTillerNumber * delt
        s.WeightLVgreen += r.dWeightLV * delt
        s.WeightLVdead += r.LVdeath * delt
        s.WeightHARV += r.dWeightHARV * delt
        s.WeightRE += r.dWeightRE * delt
        s.WeightRT += r.dWeightRT * delt
        s.LeafLength += r.dLeafLengthAct * delt

        # Total above ground dry weight including harvests, kg ha-1
        s.WeightABG = s.WeightHARV + s.WeightLVgreen

        # Running specific leaf area in model, ha kg-1
        s.SLAINT = s.LAI / s.WeightLVgreen

        s.DVS = s.TSUM / p.TSUMmax
        # TODO: decide whether to reset TSUM after mowing

        self.soil_temperature.integrate(day, delt)
        self.root_dynamics.integrate(day, delt)

    @prepare_states
    def finalize(self, day):
        SimulationObject.finalize(self, day)

    def _on_MOWING(self, biomass_remaining):
        """Handler for grass mowing events
        """
        self.WeightLV_remaining = biomass_remaining
        self._flag_MOWING = True
