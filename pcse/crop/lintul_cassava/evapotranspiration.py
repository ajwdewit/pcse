# -*- coding: utf-8 -*-
# Herman Berghuijs (herman.berghuijs@wur.nl), Allard de Wit (allard.dewit@wur.nl), Tom Schut (tom.schut@wur.nl)
# February 2026

from pcse.base import ParamTemplate, RatesTemplate, SimulationObject, StatesTemplate
from pcse.traitlets import Float
import numpy as np

class evapotranspiration(SimulationObject):
    """
    Class to simulate potential and actual transpiration and soil evaporation in LINTUL Cassava.
    
    The class first calculates the rates op potential transpiration and soil evaporation from the reference
    evapotranspiration, the rain interception and the soil moisture content. Next it calculates the
    critical soil moisture contents above which drought stress occurs or above which oxygen stress
    occurs. This is used to calculate both the actual transpiration and actual soil evaporation
    rates. Finally, a transpiration reduction factor is calculated.     
    
    **Simulation parameters**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Type     Unit
    =================  ==============================================  ======  ===========================    
    SMFCF              Soil moisture content at field capacity         SCr      cm3 water cm-2 ground
    SM0                Soil moisture content at saturation             SCr      cm3 water cm-2 ground
    SMW                Soil moisture content at wilting point          SCr      cm3 water cm-2 ground
    TRANCO             Transpiration constant that indicates the 
                       level of drought tolerance                      SCr      cm3 water cm-2 ground d-1
    TWCSD              Ratio of soil moisture content at extreme
                       drought and the soil moisture content at 
                       wilting point.                                  SCr      cm3 water cm3 water
    WCAD               Soil moisture content at air dry                SCr      cm3 water cm-3 ground
    WCWET              Soil moisture content above which oxygen
                       stress occurs                                   SCr      cm3 water cm-3 ground
    =================  ==============================================  ======  ===========================    
    
    **State variables**

    None

    **Rate variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    EVWMX              Maximum evaporation rate of an open water
                       surface                                         Y       cm3 water cm-2 ground d-1
    EVSMX              Maximum soil evaporation rate                   N       cm3 water cm-2 ground d-1
    RPEVAP             Potential soil evaporation rate                 N       cm3 water cm-2 ground d-1
    RPTRAN             Potential transpiration rate                    N       cm3 water cm-2 ground d-1
    TRA                Actual transpiration rate                       Y       cm3 water cm-2 ground d-1
    =================  ==============================================  ======  ===========================

    **Auxillary variables**

    =================  ==============================================  ======  ===========================
    Name               Description                                     Pbl     Unit
    =================  ==============================================  ======  ===========================
    RFTRA              Transpiration reduction factor                  Y       cm3 water cm-3 water
    WCCR               Critical soil moisture content below which
                       drought stress can occur                        Y       cm3 water cm-3 ground
    WCSD               Soil moisture content at severe drought         Y       cm3 water cm-3 ground
    =================  ==============================================  ======  ===========================
    """
    
    class Parameters(ParamTemplate):
        TRANCO = Float()
        TWCSD = Float()
        WCAD = Float()
        SMFCF = Float()
        SM0 = Float()
        WCWET = Float()
        SMW = Float()

    class RateVariables(RatesTemplate):
        RPEVAP = Float()
        RPTRAN = Float()
        EVSMX = Float()
        EVWMX = Float()
        TRA = Float()
        RFTRA = Float()
        WCCR = Float()
        WCSD = Float()

    class StateVariables(StatesTemplate):
        pass

    def initialize(self, day, kiosk, parameters):
        self.kiosk = kiosk
        self.params = self.Parameters(parameters)
        self.rates = self.RateVariables(kiosk,
                                        publish = ["EVSMX", "TRA", "RFTRA", "WCCR", "WCSD", "EVWMX"])
        self.states = self.StateVariables(
            kiosk,
            publish=[]
        )

    def __call__(self, day, drv, delt = 1):
        k = self.kiosk
        p = self.params
        r = self.rates

        # Potential evaporation and transpiration are weighed by a factor representing the plant canopy (exp(-0.5 * LAI)).
        RPEVAP = np.exp(-0.5 * k.LAI) * drv.ES0 # cm d-1
        RPTRAN = (1 - np.exp(-0.5 * k.LAI)) * drv.ET0 # cm d-1
        RPTRAN = max(0, RPTRAN - 0.5 * k.RNINTC)  # cm d-1

        # Evaporation is decreased when water content is below field capacity,
        # but continues until WC = WCAD. It is ensured to stay within 0-1 range
        limit_evap = (k.SM - p.WCAD) / (p.SMFCF - p.WCAD)  # (-)
        limit_evap = min(1, max(0, limit_evap))  # (-)
        EVAP = RPEVAP * limit_evap  # cm d-1

        # Maximum evaporation from an open water surface [cm d-1]. It is not used by the native soil water balance of LINTUL
        # Cassava, but other soil water balances in PCSE require this as input from evapotranspiration modules.
        REVAPW = drv.E0

        # Soil moisture content at severe drought is calculated to see if drought stress occurs in the crop.
        WCSD = p.SMW * p.TWCSD

        # Critical water content for drought stress. The critical soil moisture content depends on
        # the transpiration coefficient (TRANCO) which is a measure of how drought resistant the crop is.
        WCCR = p.SMW + max(WCSD - p.SMW, RPTRAN / (RPTRAN + p.TRANCO) * (p.SMFCF - p.SMW))

        # Reduction factor RF for transpiration can be either from drought or oxygen stress
        if k.SM <= p.SMW:
            # soil moisture below wilting point: no transpiration possible
            FR = 0.0
        elif p.SMW < k.SM <= WCCR:
            # Soil moisture between wilting point and critical level: reduced transpiration
            FR = (k.SM - p.SMW) / (WCCR - p.SMW)
        elif WCCR < k.SM <= p.WCWET:
            # soil moisture above critical level and below wet level: no reduction
            FR = 1.0
        else:
            # soil moisture above wet level: reduction due to oxygen stress
            FR = (p.SM0 - k.SM) / (p.SM0 - p.WCWET)

        # Ensure to stay within the 0-1 range
        FR = min(1.0, max(0.0, FR))

        # Actual transpiration
        TRAN = RPTRAN * FR  # cm d-1

        # A final correction term is calculated to reduce evaporation and transpiration when
        # evapotranspiration exceeds the amount of water in soil present in excess of air dryness.
        WAAD = p.WCAD * k.RD  # The amount of soil water at air dryness (AD) [cm]
        W_avail = k.SM * k.RD - WAAD # Actual available amount in access of AD [cm]
        W_required = (EVAP + TRAN) * delt  # Amount of evapotranspiration [cm]
        if W_required > W_avail:  # more water is asked than available in the soil -> reduce ET rates.
            AVAILF = W_avail / W_required
            TRAN = TRAN * AVAILF
            EVAP = EVAP * AVAILF

        # The final transpiration reduction factor is defined as the ratio between actual and
        # potential transpiration
        RFTRA = TRAN / RPTRAN if RPTRAN > 0 else 1.0

        r.RPEVAP = RPEVAP
        r.RPTRAN = RPTRAN

        r.EVWMX = REVAPW
        r.EVSMX = EVAP
        r.TRA = TRAN
        r.RFTRA = RFTRA
        r.WCCR = WCCR
        r.WCSD = WCSD