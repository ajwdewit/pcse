# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl) and Herman Berghuijs (herman.berghuijs@wur.nl), January 2024
"""This module wraps the soil components for water and nutrients so that they run jointly
within the same model. Use of these wrappers is only required for WOFOST 8+ because the
WOFOST 7x version do not simulated nutrient-limited production and therefore can import
the waterbalance directly into the configuration.
"""
from pcse.base import SimulationObject
from .classic_waterbalance import WaterbalanceFD, WaterbalancePP
from .multilayer_waterbalance import WaterBalanceLayered, WaterBalanceLayered_PP
from .n_soil_dynamics import N_Soil_Dynamics, N_PotentialProduction
from .snomin import SNOMIN
from ..traitlets import Instance


class BaseSoilWrapper(SimulationObject):
    """Base class for wrapping soil water and nutrient/carbon balances.
    """
    waterbalance_class = None
    nutrientbalance_class = None
    waterbalance = Instance(SimulationObject)
    nutrientbalance = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        if self.waterbalance_class is not None:
            self.waterbalance = self.waterbalance_class(day, kiosk, parvalues)
        if self.nutrientbalance_class is not None:
            self.nutrientbalance = self.nutrientbalance_class(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        if self.waterbalance_class is not None:
            self.waterbalance.calc_rates(day, drv)
        if self.nutrientbalance_class is not None:
            self.nutrientbalance.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        if self.waterbalance_class is not None:
            self.waterbalance.integrate(day, delt)
        if self.nutrientbalance_class is not None:
            self.nutrientbalance.integrate(day, delt)


class SoilModuleWrapper_PP(BaseSoilWrapper):
    """This wraps the soil water balance and soil N balance for potential production.
    """
    waterbalance_class = WaterbalancePP
    nutrientbalance_class = N_PotentialProduction


class SoilModuleWrapper_WLP_CWB(BaseSoilWrapper):
    """This wraps the classic soil water balance for free drainage conditions and unlimited N balance
    for production conditions limited by soil water only.
    """
    waterbalance_class = WaterbalanceFD
    nutrientbalance_class = N_PotentialProduction


class SoilModuleWrapper_NWLP_CWB_CNB(BaseSoilWrapper):
    """This wraps the classic soil water balance and classic N balance
    for production conditions limited by both soil water and N but with
    simple water and N dynamics.
    """
    waterbalance_class = WaterbalanceFD
    nutrientbalance_class = N_Soil_Dynamics


class SoilModuleWrapper_WLP_MLWB(BaseSoilWrapper):
    """This wraps the multi-layer soil water balance and classic N balance
    for production conditions limited by both soil water and N.
    """
    waterbalance_class = WaterBalanceLayered
    nutrientbalance_class = N_PotentialProduction

class SoilModuleWrapper_NWLP_MLWB_CNB(BaseSoilWrapper):
    """This wraps the soil water balance for free drainage conditions and N balance
    for production conditions limited by both soil water and N but with
    advanced soil water dynamics and simple N dynamics.
    """
    waterbalance_class = WaterBalanceLayered
    nutrientbalance_class = N_Soil_Dynamics


class SoilModuleWrapper_NWLP_MLWB_SNOMIN(BaseSoilWrapper):
    """This wraps the soil water balance for free drainage conditions and the
    advanced SNOMIN C/N balance for production conditions limited by both soil water
    and N using advanced soil water dynamics and C/N dynamics.
    """
    waterbalance_class = WaterBalanceLayered
    nutrientbalance_class = SNOMIN

