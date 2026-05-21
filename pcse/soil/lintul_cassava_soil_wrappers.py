# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl) and Herman Berghuijs (herman.berghuijs@wur.nl), January 2024
"""This module wraps the soil components for water and nutrients so that they run jointly
within the same model.
"""
from pcse.base import SimulationObject
from pcse.traitlets import Instance
from pcse.soil.classic_waterbalance import WaterbalanceFD, WaterbalancePP
from pcse.soil.lintul_cassava_soil_nutrient_dynamics import soil_nutrient_dynamics, soil_nutrient_dynamics_PP


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


class Lintul_Cassava_NWLP_SoilWrapper(BaseSoilWrapper):
    waterbalance_class = WaterbalanceFD
    nutrientbalance_class = soil_nutrient_dynamics


class Lintul_Cassava_WLP_SoilWrapper(BaseSoilWrapper):
    waterbalance_class = WaterbalanceFD
    nutrientbalance_class = soil_nutrient_dynamics_PP


class Lintul_Cassava_NLP_SoilWrapper(BaseSoilWrapper):
    waterbalance_class = WaterbalancePP
    nutrientbalance_class = soil_nutrient_dynamics


class Lintul_Cassava_PP_SoilWrapper(BaseSoilWrapper):
    waterbalance_class = WaterbalancePP
    nutrientbalance_class = soil_nutrient_dynamics_PP
