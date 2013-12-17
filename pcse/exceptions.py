"""Exception hierarchy for PCSE
"""

class PCSEError(Exception):
    """Top PyWOFOST Exception"""
    
class CarbonBalanceError(PCSEError):
    "Raised when carbon flows are not balanced."

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