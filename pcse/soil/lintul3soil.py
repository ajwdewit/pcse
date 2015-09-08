# -*- coding: utf-8 -*-
from ..base_classes import ParamTemplate, SimulationObject
from ..base_classes import StatesWithImplicitRatesTemplate as StateVariables
from ..decorators import prepare_rates, prepare_states
from ..traitlets import Float, Bool
from ..util import limit
from math import sqrt
from ..exceptions import WaterBalanceError

cm2mm = lambda cm: 10. * cm
m2mm = lambda x: 1000 * x

class Lintul3Soil(SimulationObject):
    """
        * ORIGINAL COPYRIGHT NOTICE:
        *-------------------------------------------------------------------------*
        * Copyright 2013. Wageningen University, Plant Production Systems group,  *
        * P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
        * You may not use this work except in compliance with the Licence.        *
        * You may obtain a copy of the Licence at:                                *
        *                                                                         *
        * http://models.pps.wur.nl/content/licence-agreement                      *
        *                                                                         *
        * Unless required by applicable law or agreed to in writing, software     *
        * distributed under the Licence is distributed on an "AS IS" basis,       *
        * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
        *-------------------------------------------------------------------------*
        
        The water balance in the model does not deal with a flooded rice system, but the 
        soil is maintained at saturation so that the crop does not experience water 
        stress. The soil water balance is calculated for a single soil layer, whose 
        thickness increases with downward root elongation. The model does not simulate 
        root elongation, but only average root increase in depth. The assumption that 
        the root elongation does not depend on the dryness of the soil until wilting 
        point is originally from LINTUL1 (Spitters, 1990). In the rice situation, we 
        have a saturated soil and the wilting point is not expected occur. Water and 
        nutrient uptake by the crop is limited to the rooted soil depth. Addition of 
        water to the available store through root extension is calculated from the rate 
        of extension and the (saturated) water content (Spitters and Schapendonk, 1990; 
        Farr� et al., 2000). 
        
        Nitrogen supply from soil 
        
        Mineral nitrogen available for crop uptake (TNSOIL) originates from three 
        sources: nitrogen present in the soil profile at germination/transplanting, 
        nitrogen from biological fixation and mineralization from soil organic matter 
        during the growing season and nitrogen applied as fertilizer. Under aerobic 
        conditions, indigenous soil nitrogen supply can be quantified reasonably 
        accurately on the basis of soil organic matter content (Sinclair and Amir, 
        1992). Under anaerobic conditions, however, differences in mineral N supply 
        between fields or seasons could not be explained on the basis of soil organic 
        carbon, total nitrogen or initial inorganic nitrogen (Cassman et al., 1996; 
        Bouman et al., 2001). Hence, indigenous nitrogen supply is introduced as a 
        site-specific exogenous input, rather than simulating nitrogen mineralization. 
        Ten Berge et al. (1997) found indigenous nitrogen supply values for tropical 
        soils in the range of 0.5�0.9 kg ha-1 d-1. Fertilizer nitrogen available for 
        plant uptake, taking into account possible losses due to volatilization, 
        denitrification, and leaching, is included in the model as a variable fraction, 
        named nitrogen recovery fraction, NRF. A low N recovery value in flooded rice 
        systems (30�39%, Cassman et al., 1996) is mainly due to volatilization losses. 
        Any form of inorganic nitrogen (NH4 + or NO3 -) when available in free form, if 
        not taken up by the crop, is vulnerable to loss either through volatilization, 
        denitrification or leaching. On a daily basis, the average N recovery value 
        would be more or less stable under flooded conditions. The variable NRF used in 
        the model represents the net recovery fraction of applied fertilizer, which 
        depends on soil type, development stage of the crop, fertilizer type and time 
        and mode of nitrogen application (De Datta, 1986) is determined by calibration.


        See Farré, I., Van Oijen, M., Leffelaar, P.A., Faci, J.M., 2000. Analysis of maize growth
            for different irrigation strategies in northeastern Spain. Eur. J. Agron. 12, 225–238.
            doi:10.1016/S1161-0301(00)00051-4
            
            
        *parameters*
        ======== =============================================== =======  ==========
         Name     Description                                     Type     Unit
        ======== =============================================== =======  ==========
        DRATE    Maximum drainage rate of the soil                          mm/day
        IRRIGF   Irrigation switch                                          (bool)
        WCAD     Water content at air dryness                               m³/m³
        WCFC     Water content at field capacity (0.03 MPa)                 m³/m³
        WCST     Water content at full saturation                           m³/m³
        WCSUBS   water content subsoil (?)                                  m³/m³
        WCWP     Water content at wilting point (1.5 MPa)                   m³/m³
        WMFAC    water management (False = irrigated up to the 
                 field capacity, true= irrigated up to saturation)          (bool)  
        ROOTDI   Initial rooting depth                                      m
        WCI      Initial water content in soil                              m³/³
        
        **State variables:**
        =========== ================================================= ==== ===============
         Name        Description                                      Pbl      Unit
        =========== ================================================= ==== ===============
        WA          Soil water amount                                  *     mm
        WC          Volumetric soil water content                            mm
        TRUNOF      Run off accumulator                                      mm
        TTRAN       Crop transpiration accumulated over growth period        mm
        TEVAP       Soil evaporation accumulated over growth period          mm
        TDRAIN      Drained accumulator                                      mm
        TRAIN       Rain accumulator                                         mm
        TEXPLO      Exploration accumulator                                  mm
        TIRRIG      Irrigation accumulator                                   mm        
    """

    class Parameters(ParamTemplate):
        DRATE   = Float(-99)    # Maximum drainage rate of the soil (mm/day)
        IRRIGF  = Bool()        # Irrigation factor 
        WCFC    = Float(-99)    # Soil hydraulic properties
        WCI     = Float(-99)    # Initial water content in cm3 of water/(cm3 of soil).
        WCST    = Float(-99)    # Soil hydraulic properties
        WCSUBS  = Float(-99)    # water content subsoil (?)
        WCAD    = Float(-99)    # Water content at air dryness m3/ m3
        WCWP    = Float(-99)    # Soil hydraulic properties
        WMFAC   = Bool()        # water management (0=irrigated up to the field capacity, 1 = irrigated up to saturation)        
        ROOTDI  = Float(-99)    # initial rooting depth [m] 

    class Lintul3SoilStates(StateVariables):
        WA      = Float(-99.)   # Amount of soil water
        WC      = Float(-99.)   # Volumetric soil water content
        TRUNOF  = Float(-99.)   # total run off
        TTRAN   = Float(-99.)   # Crop transpiration accumulated over growth period
        TEVAP   = Float(-99.)   # soil evaporation accumulated over growth period
        TDRAIN  = Float(-99.)   # todtal drained
        TRAIN   = Float(-99.)   # total rain
        TEXPLO  = Float(-99.)   # Total Exploration
        TIRRIG  = Float(-99.)   # total irrigation

    # Counter for days since last rain
    DSLR  = 0
    # Initial amount of water in upper layer (rooted zone)
    WAI = 0.

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)

        init = self.Lintul3SoilStates.initialValues()

        # Initial amount of water present in the rooted depth at the start of
        # the calculations, based on the initial water content (in mm).
        self.WAI = m2mm(self.params.ROOTDI) * self.params.WCI
        init["WA"] = self.WAI
        
        # Initialize state variables
        self.states = self.Lintul3SoilStates(kiosk, publish=["WA", "WC"], **init)
        self.states.initialize_rates()

    def _safe_get_from_kiosk(self, varname, default=0.0):
        """
        Get named value from the kiosk; return default if it isn't available
        :param varname:    variable name
        :param default:    returned if the value is not available from the kiosk
                           (optional, default 0.) 
        """
        if varname in self.kiosk:
            return self.kiosk[varname]
        return default

    @prepare_rates
    def calc_rates(self, day, drv):

        # dynamic calculations
        p = self.params
        s = self.states

        DELT = 1.  #

        ROOTD = self._safe_get_from_kiosk("ROOTD", p.ROOTDI)
        RROOTD = self._safe_get_from_kiosk("RROOTD")
        PEVAP = self._safe_get_from_kiosk("PEVAP", cm2mm(drv.ES0))
        TRAN = self._safe_get_from_kiosk("TRAN")
        
        # Variables supplied by the weather system
        RAIN = cm2mm(drv.RAIN)  # cm  --> mm CORRECTION FOR NON-STANDARD cm in WOFOST-WEATHER

        # Soil evaporation rate
        EVAP = self._soil_evaporation(RAIN, PEVAP, ROOTD, DELT)
        
        # Calling the subroutine for rates of drainage, runoff and irrigation.
        DRAIN, RUNOFF, IRRIG = self._drainage_runoff_irrigation(RAIN, EVAP, TRAN, DELT, s.WA, ROOTD)

        #  Exploration of water in soil when roots grow downward.
        EXPLOR = m2mm(RROOTD) * p.WCSUBS
                
        RWA = (RAIN + EXPLOR + IRRIG)-(RUNOFF + TRAN + EVAP + DRAIN)

        # Assign rate variables for associated states
        s.rWA = RWA
        s.rTEXPLO = EXPLOR
        s.rTEVAP = EVAP
        s.rTTRAN = TRAN
        s.rTRUNOF = RUNOFF
        s.rTIRRIG = IRRIG 
        s.rTRAIN  = RAIN  
        s.rTDRAIN = DRAIN
        
        WATBAL = (s.WA + (s.TRUNOF + s.TTRAN + s.TEVAP + s.TDRAIN)
                       - (self.WAI + s.TRAIN + s.TEXPLO + s.TIRRIG))

        if abs(WATBAL) > 0.0001:
            raise WaterBalanceError("water un-balance in root zone at day %s" % day)

    @prepare_states
    def integrate(self, day):
        s = self.states
        p = self.params
        s.integrate(delta=1.)

        # Volumetric Water content in the rootzone
        ROOTD = self._safe_get_from_kiosk("ROOTD", p.ROOTDI)
        s.WC = s.WA / m2mm(ROOTD)

    def _drainage_runoff_irrigation(self, RAIN, EVAP, TRAN, DELT, WA, ROOTD):
        """compute rates of drainage, runoff and irrigation
        """
        p = self.params

        WAFC = p.WCFC * m2mm(ROOTD)
        WAST = p.WCST * m2mm(ROOTD)

        DRAIN = limit(0., p.DRATE, (WA-WAFC)/DELT + (RAIN - EVAP - TRAN))
        RUNOFF = max(0., (WA-WAST)/DELT + (RAIN - EVAP - TRAN - DRAIN))
        
        if p.WMFAC:
            # If a soil is irrigated by flooding, : soil water content is
            # kept at saturation via "irrigation events".
            IRRIG = max(0., (WAST-WA)/DELT - (RAIN - EVAP - TRAN - DRAIN - RUNOFF)) if p.IRRIGF else 0.0
        else:
            # If soil is irrigated but not flooded, : soil water content
            # is kept at field capacity via "irrigation events".
            IRRIG = max(0., (WAFC-WA)/DELT - (RAIN - EVAP - TRAN - DRAIN - RUNOFF)) if p.IRRIGF else 0.0
    
        return DRAIN, RUNOFF, IRRIG

    def _soil_evaporation(self, RAIN, PEVAP, ROOTD, DELT):
        """Compute actual soil evaporation rate as function
        of Days-Since-Last-Rain (DSLR) limiting for the amount
        of water at air-dry [WAAD]

        :param RAIN: Rainfall amount in [mm]
        :param PEVAP: Potential soil evaporation rate [mm/d]
        :param DELT: time step [1 day]
        """
        
        p = self.params
        s = self.states
        
        # see also classic_water balance.py
        WAAD = p.WCAD * m2mm(ROOTD)
         
        if RAIN >= 0.5:
            EVS = PEVAP
            self.DSLR = 1.
        else:
            self.DSLR += 1.
            EVSMXT = PEVAP*(sqrt(self.DSLR) - sqrt(self.DSLR - 1.))
            EVS = min(PEVAP, EVSMXT + RAIN)
        
        # WA-WAAD is the physically available amount of water in the soil above
        # air-dry. Limit the amount of evapotranspiration for this amount  
        # in order to avoid emptying the soil water reservoir below air-dry
        AVAILF = min(1., (s.WA - WAAD)/(EVS * DELT)) if (EVS > 0) else 0.0
        
        return EVS * AVAILF
