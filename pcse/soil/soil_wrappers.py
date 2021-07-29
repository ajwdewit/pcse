# -*- coding: utf-8 -*-
# Copyright (c) 2004-2020 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), September 2020
"""This module wraps the soil components for water and nutrients so that they run jointly
within the same model.
"""
from pcse.base import SimulationObject
from .classic_waterbalance import WaterbalanceFD, WaterbalancePP
from .npk_soil_dynamics import NPK_Soil_Dynamics, NPK_PotentialProduction
from .n_soil_dynamics import N_PotentialProduction, N_Soil_Dynamics
from ..traitlets import Instance
from ..decorators import prepare_states


class SoilModuleWrapper_PP(SimulationObject):
    """This wraps the soil water balance and soil NPK balance for potential production.
    """
    WaterbalancePP = Instance(SimulationObject)
    NPK_PotentialProduction = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        self.WaterbalancePP = WaterbalancePP(day, kiosk, parvalues)
        self.NPK_PotentialProduction = NPK_PotentialProduction(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.WaterbalancePP.calc_rates(day, drv)
        self.NPK_PotentialProduction.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        self.WaterbalancePP.integrate(day, delt)
        self.NPK_PotentialProduction.integrate(day, delt)


class SoilModuleWrapper_WLP_FD(SimulationObject):
    """This wraps the soil water balance for free drainage conditions and NPK balance
    for production conditions limited by soil water only.
    """
    WaterbalanceFD = Instance(SimulationObject)
    NPK_Soil_Dynamics = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        self.WaterbalanceFD = WaterbalanceFD(day, kiosk, parvalues)
        self.NPK_Soil_Dynamics = NPK_PotentialProduction(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.WaterbalanceFD.calc_rates(day, drv)
        self.NPK_Soil_Dynamics.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        self.WaterbalanceFD.integrate(day, delt)
        self.NPK_Soil_Dynamics.integrate(day, delt)


class SoilModuleWrapper_NPK_WLP_FD(SimulationObject):
    """This wraps the soil water balance for free drainage conditions and NPK balance
    for production conditions limited by both soil water and NPK.
    """
    WaterbalanceFD = Instance(SimulationObject)
    NPK_Soil_Dynamics = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        self.WaterbalanceFD = WaterbalanceFD(day, kiosk, parvalues)
        self.NPK_Soil_Dynamics = NPK_Soil_Dynamics(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.WaterbalanceFD.calc_rates(day, drv)
        self.NPK_Soil_Dynamics.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        self.WaterbalanceFD.integrate(day, delt)
        self.NPK_Soil_Dynamics.integrate(day, delt)


class SoilModuleWrapper_N_WLP_FD(SimulationObject):
    """This wraps the soil water balance for free drainage conditions and N balance
    for production conditions limited by both soil water and N.
    """
    WaterbalanceFD = Instance(SimulationObject)
    N_Soil_Dynamics = Instance(SimulationObject)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with parameter key/value pairs
        """
        self.WaterbalanceFD = WaterbalanceFD(day, kiosk, parvalues)
        self.N_Soil_Dynamics = N_Soil_Dynamics(day, kiosk, parvalues)

    def calc_rates(self, day, drv):
        self.WaterbalanceFD.calc_rates(day, drv)
        self.N_Soil_Dynamics.calc_rates(day, drv)

    def integrate(self, day, delt=1.0):
        self.WaterbalanceFD.integrate(day, delt)
        self.N_Soil_Dynamics.integrate(day, delt)
