# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl) and Herman Berghuijs (herman.berghuijs@wur.nl), April 2024
"""Exception hierarchy for PCSE
"""

class PCSEError(Exception):
    """Top PCSE Exception"""
    
class CarbonBalanceError(PCSEError):
    "Raised when carbon flows are not balanced."

class NutrientBalanceError(PCSEError):
    "Raised when nitrogen flows are not balanced."

class WaterBalanceError(PCSEError):
    "Raised when water balance does not close"
    
class PartitioningError(PCSEError):
    "Raised when problems with partitioning are found."
    
class ParameterError(PCSEError):
    "Raised when problems with parameters are found."
    
class VariableKioskError(PCSEError):
    "Raised when problems with kiosk registrations are found."
    
class WeatherDataProviderError(PCSEError):
    "Raised when problems occur with the WeatherDataProviders" 

class SoilOrganicMatterBalanceError(PCSEError):
    "Raised when soil organic matter balance does not close"

class SoilOrganicCarbonBalanceError(PCSEError):
    "Raised when soil organic carbon balance does not close"

class SoilOrganicNitrogenBalanceError(PCSEError):
    "Raised when soil organic nitrogen balance does not close"

class SoilAmmoniumBalanceError(PCSEError):
    "Raised when soil ammonium balance does not close"

class SoilNitrateBalanceError(PCSEError):
    "Raised when soil nitrate balance does not close"