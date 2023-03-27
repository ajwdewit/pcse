# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
import numpy as np
from .. import exceptions as exc
from pcse.traitlets import Float
from pcse.decorators import prepare_rates, prepare_states
from pcse.base import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from pcse import signals
from ..traitlets import Float, Int, Instance, Bool
from .soiln_profile import SoilNProfile

class N_soil_dynamics_layered(SimulationObject):
    """Provides unlimited soil N/P/K for potential production simulations.

    NAVAIL just remains 100 kg/ha whatever the crop takes.
    """

    # Placeholders initial values
    _ORGMATI = None
    _CORGI = None
    _NORGI = None
    _NH4I = None
    _NO3I = None

    # Placeholders
    _RNO3AM = None
    _RNH4AM = None
    _RAGEAM = None
    _RORGMATAM = None
    _RCORGAM = None
    _RNORGAM = None

    # Unit conversions
    g_to_kg = 1e-3
    cm_to_m = 1e-2
    cm2_to_ha = 1e-8
    cm3_to_m3 = 1e-6
    ha_to_m2 = 1e-4
    m2_to_ha = 1e-4
    y_to_d = 365.25

    # placeholder for soil object
    soiln_profile = None

    class StateVariables(StatesTemplate):
        AGE0   = Instance(np.ndarray) # Initial age of material (d)
        AGE    = Instance(np.ndarray) # Appearant age of material (d)
        ORGMAT = Instance(np.ndarray) # Amount of organic matter (kg ORG ha-1)
        CORG   = Instance(np.ndarray) # Amount of C in organic matter (kg C ha-1)
        NORG   = Instance(np.ndarray) # 
        NH4    = Instance(np.ndarray) # Amount of NH4-N (kg N ha-1)
        NO3    = Instance(np.ndarray) # Amount of NO3-N (kg N ha-1)
        NAVAIL = Float()  # total mineral N from soil and fertiliser  kg N ha-1
        NDENITCUM = Float()
        NO3LEACHCUM = Float()
        NH4LEACHCUM = Float()
        NLOSSCUM = Float()

        RORGMATDISTT = Float()
        RORGMATAMTT = Float()
        RCORGDISTT = Float()
        RCORGAMTT = Float()
        RNORGDISTT = Float()
        RNORGAMTT = Float()

        RNO3NITRTT = Float()
        RNO3DENITRTT = Float()
        RNO3UPTT = Float()
        RNO3INTT = Float()
        RNO3OUTTT = Float()
        RNO3AMTT = Float()

        RNH4MINTT = Float()
        RNH4NITRTT = Float()
        RNH4UPTT = Float()
        RNH4INTT = Float()
        RNH4OUTTT = Float()
        RNH4AMTT = Float()

        ORGMATT = Float()
        CORGT = Float()       
        NORGT = Float()
        RMINT = Float()
        NH4T = Float()
        NO3T = Float()

    class RateVariables(RatesTemplate):
        RAGE = Instance(np.ndarray)
        RORGMAT = Instance(np.ndarray)
        RCORG = Instance(np.ndarray)
        RNORG = Instance(np.ndarray)

        RAGEAG = Instance(np.ndarray)

        RORGMATDIS = Instance(np.ndarray)
        RCORGDIS = Instance(np.ndarray)
        RNORGDIS = Instance(np.ndarray)

        RAGEAM = Instance(np.ndarray)
        RORGMATAM = Instance(np.ndarray)
        RCORGAM = Instance(np.ndarray)
        RNORGAM = Instance(np.ndarray)
        
        RNH4 = Instance(np.ndarray)
        RNH4MIN = Instance(np.ndarray)
        RNH4NITR = Instance(np.ndarray)
        RNH4UP = Instance(np.ndarray)
        RNH4IN = Instance(np.ndarray)
        RNH4OUT = Instance(np.ndarray)
        RNH4AM = Instance(np.ndarray)

        RNO3 = Instance(np.ndarray)
        RNO3NITR = Instance(np.ndarray)
        RNO3DENITR = Instance(np.ndarray)
        RNO3UP = Instance(np.ndarray)
        RNO3IN = Instance(np.ndarray)
        RNO3OUT = Instance(np.ndarray)
        RNO3AM = Instance(np.ndarray)

        RNH4LEACHCUM = Float()
        RNO3LEACHCUM = Float()
        RNDENITCUM = Float()
        RNLOSS = Float()

    class Parameters(ParamTemplate):
        A0SOM = Float()             # Initial age of humus (y)
        CNRatioBio = Float()        # C:N Ratio of microbial biomass (kg C kg-1 N)
        FASDIS = Float()            # Fraction of assimilation to dissimilation (kg ORG kg-1 ORG)
        KDENIT_REF = Float()        # Reference first order denitrification rate constant (d-1)
        KNIT_REF = Float()          # Reference first order nitrification rate constant (d-1)
        KSORP = Float()             # Sorption coefficient ammonium (m3 water kg-1 soil)
        MRCDIS = Float()            # Michaelis Menten constant for response factor denitrification to soil respiration
        NO3I = Instance(list)       # Initial amount of NO3-N (kg N ha-1)
        NH4I = Instance(list)       # Initial amount of NH4-N (kg N ha-1)
        WFPS_CRIT = Float()         # Critical water filled pore space fraction (m3 water m-3 pore) for denitrification

    def initialize(self, day, kiosk, parvalues):
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)
        self.soiln_profile = SoilNProfile(parvalues)

        # Initialize module
        sinm = self.SoilInorganicNModel()

        # Initialize state variables
        parvalues._soildata["soil_profile"] = self.soiln_profile
        NH4 = np.zeros(len(self.soiln_profile))
        NO3 = np.zeros_like(NH4)
        AGE =  np.zeros((1, len(self.soiln_profile)))
        AGE0 = np.zeros_like(AGE)
        ORGMAT = np.zeros_like(AGE)
        CORG =  np.zeros_like(AGE)
        NORG =  np.zeros_like(AGE)
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            NH4[il] = self.params.NH4I[il] * self.m2_to_ha
            NO3[il] = self.params.NO3I[il] * self.m2_to_ha
            AGE0[0,il] = self.params.A0SOM * self.y_to_d
            AGE[0,il] = self.params.A0SOM * self.y_to_d
            ORGMAT[0,il] = layer.RHOD_kg_per_m3 * layer.FSOMI * layer.Thickness_m
            CORG[0,il] = minip_C.calculate_organic_C(ORGMAT[0,il])
            NORG[0,il] = CORG[0, il] / layer.CNRatioSOMI

        # Initialize state variables to check the mass balance of organic matter
        RORGMATDISTT = 0.
        RORGMATAMTT = 0.

        # Initialize state variables to check the mass balance of organic C
        RCORGDISTT = 0.
        RCORGAMTT = 0.

        # Initialize state variables to check the mass balance of organic N
        RNORGDISTT = 0.
        RNORGAMTT = 0.

        # Initialize state variables to check the mass balance of NH4-N
        RNH4MINTT = 0.
        RNH4NITRTT = 0.
        RNH4UPTT = 0.
        RNH4INTT = 0.
        RNH4OUTTT = 0.
        RNH4AMTT = 0.

        # Initialize state variables to check the mass balance of NO3-N
        RNO3NITRTT = 0.
        RNO3DENITRTT = 0.
        RNO3UPTT = 0.
        RNO3INTT = 0.
        RNO3OUTTT = 0.
        RNO3AMTT = 0.

        # Initialize output variables
        CORGT = np.sum(CORG) / self.m2_to_ha
        NORGT = np.sum(NORG) / self.m2_to_ha
        ORGMATT = np.sum(ORGMAT) / self.m2_to_ha
        RMINT = 0.
        NH4T = np.sum(NH4)/ self.m2_to_ha
        NO3T = np.sum(NO3)/ self.m2_to_ha
        NAVAIL = 0.
        NH4LEACHCUM = 0.
        NO3LEACHCUM = 0.
        NDENITCUM = 0.
        NLOSSCUM = 0.

        states = {"NAVAIL": NAVAIL, "NO3": NO3, "NH4": NH4, "AGE": AGE, "AGE0": AGE0, "ORGMAT": ORGMAT, "CORG": CORG, "NORG": NORG,  
                  "ORGMATT": ORGMATT,  "CORGT": CORGT, "NORGT": NORGT, "RMINT": RMINT, "NH4T": NH4T, "NO3T": NO3T, 
                  "RORGMATAMTT": RORGMATAMTT, "RORGMATDISTT": RORGMATDISTT, "RCORGAMTT": RCORGAMTT, "RCORGDISTT": RCORGDISTT,
                  "RNORGAMTT": RNORGAMTT, "RNORGDISTT": RNORGDISTT, "RNH4MINTT": RNH4MINTT, "RNH4NITRTT": RNH4NITRTT, "RNH4UPTT": RNH4UPTT,
                  "RNH4INTT": RNH4INTT, "RNH4OUTTT": RNH4OUTTT, "RNH4AMTT": RNH4AMTT, "RNO3NITRTT": RNO3NITRTT, "RNO3DENITRTT": RNO3DENITRTT,
                  "RNO3UPTT": RNO3UPTT, "RNO3INTT": RNO3INTT, "RNO3OUTTT": RNO3OUTTT, "RNO3AMTT": RNO3AMTT,
                  "NH4LEACHCUM": NH4LEACHCUM, "NO3LEACHCUM": NO3LEACHCUM, "NDENITCUM": NDENITCUM, "NLOSSCUM": NLOSSCUM}

        self.states = self.StateVariables(kiosk, publish=["NAVAIL", "ORGMATT", "CORGT", "NORGT"], **states)
        self.rates = self.RateVariables(kiosk)

        k = self.kiosk
        s = self.states
        p = self.params

        # Get external state variable (in meters)
        RD_m = self.get_root_length(k)
        SM = self.get_soil_moisture_content(k)        
        NAVAIL = sinm.calculate_NAVAIL(self.soiln_profile, p.KSORP, s.NH4, s.NO3, RD_m, SM) / self.m2_to_ha

        self._RAGEAM = np.zeros_like(AGE)
        self._RORGMATAM = np.zeros_like(ORGMAT)
        self._RCORGAM = np.zeros_like(CORG)
        self._RNORGAM = np.zeros_like(NORG)
        self._RNH4AM = np.zeros_like(NH4)
        self._RNO3AM = np.zeros_like(NO3)

        self._ORGMATI = ORGMAT
        self._CORGI = CORG
        self._NORGI = NORG
        self._NH4I = NH4
        self._NO3I = NO3

        # Connect module to signal AgroManager
        self._connect_signal(self._on_APPLY_N, signals.apply_n)

    @prepare_rates
    def calc_rates(self, day, drv):
        k = self.kiosk
        r = self.rates
        s = self.states
        p = self.params

        delt = 1.0

        # Initialize model components
        sonm = self.SoilOrganicNModel()
        sinm = self.SoilInorganicNModel()

        # Get external variables in the right unit
        flow_m_per_d = self.get_water_flow_rates(k)
        N_demand_soil = self.get_N_demand(k)
        RD_m = self.get_root_length(k)
        SM = self.get_soil_moisture_content(k)
        T = drv.TEMP

        # Collect ammendment rates
        r.RAGEAM = self._RAGEAM
        r.RORGMATAM = self._RORGMATAM
        r.RCORGAM = self._RCORGAM
        r.RNORGAM = self._RNORGAM
        r.RNH4AM = self._RNH4AM
        r.RNO3AM = self._RNO3AM

        # Reset placeholders for ammendment rates
        self._RAGEAM = np.zeros_like(s.AGE)
        self._RORGMATAM = np.zeros_like(r.RORGMATAM)
        self._RCORGAM = np.zeros_like(r.RCORGAM)
        self._RNORGAM = np.zeros_like(r.RNORGAM)
        self._RNH4AM = np.zeros_like(r.RNH4)
        self._RNO3AM = np.zeros_like(r.RNH4)

        # Calculate increase aparent age of each ammendment
        r.RAGEAG = sonm.calculate_apparent_age_increase_rate(s.AGE, delt, T)

        # Calculate dissimilation rates of each ammendment
        r.RORGMATDIS, r.RCORGDIS, r.RNORGDIS = sonm.calculate_dissimilation_rates(s.AGE, s.AGE0, p.CNRatioBio, p.FASDIS, s.NORG, s.ORGMAT, T)

        # Calculate rates of change apparent age, organic matter, organic C, and organic N
        r.RAGE = r.RAGEAG + r.RAGEAM
        r.RORGMAT = r.RORGMATAM - r.RORGMATDIS
        r.RCORG = r.RCORGAM - r.RCORGDIS
        r.RNORG = r.RNORGAM - r.RNORGDIS

        # Calculate N uptake rates
        r.RNH4UP, r.RNO3UP = sinm.calculate_N_uptake_rates(self.soiln_profile, delt, p.KSORP, N_demand_soil, s.NH4, s.NO3, RD_m, SM)

        # Calculate remaining amounts of NH4-N and NO3-N after uptake and calculate chemical conversion
        NH4PRE = s.NH4 - r.RNH4UP.sum() * delt
        NO3PRE = s.NO3 - r.RNO3UP.sum() * delt
        r.RNH4MIN, r.RNH4NITR, r.RNO3NITR, r.RNO3DENITR =  sinm.calculate_reaction_rates(self.soiln_profile, p.KDENIT_REF, p.KNIT_REF, 
                                                                                         p.KSORP, p.MRCDIS, NH4PRE, NO3PRE, r.RCORGDIS, r.RNORGDIS, SM, T, p.WFPS_CRIT)
        
        # Calculate remaining amounts of NH4-N and NO3-N after uptake and reaction and calculate inorganic N flow rates between layers
        NH4PRE2 = NH4PRE + (r.RNH4MIN.sum() - r.RNH4NITR.sum()) * delt
        NO3PRE2 = NO3PRE + (r.RNO3NITR.sum() - r.RNO3DENITR.sum()) * delt
        r.RNH4IN, r.RNH4OUT, r.RNO3IN, r.RNO3OUT = sinm.calculate_flow_rates(self.soiln_profile, flow_m_per_d, p.KSORP, NH4PRE2, NO3PRE2, SM)

        # Calculate rates of change NH4-N and NO3-N
        r.RNH4 = r.RNH4AM + r.RNH4MIN - r.RNH4NITR - r.RNH4UP + r.RNH4IN - r.RNH4OUT
        r.RNO3 = r.RNO3AM + r.RNO3NITR - r.RNO3DENITR - r.RNO3UP + r.RNO3IN - r.RNO3OUT

        # Calculate output rate variables for N loss
        r.RNH4LEACHCUM =  (1/self.m2_to_ha) * r.RNH4OUT[-1]
        r.RNO3LEACHCUM =  (1/self.m2_to_ha) * r.RNO3OUT[-1]
        r.RDENITCUM = (1/self.m2_to_ha) * r.RNO3DENITR.sum()
 
    @prepare_states
    def integrate(self, day, delt=1.0):
        k = self.kiosk
        p = self.params
        r = self.rates
        s = self.states

        # Initialize soil module
        sinm = self.SoilInorganicNModel()

        # Calculate apparent ages and amounts of organic matter, organic C, and organic N in next time step
        AGE = s.AGE + r.RAGE * delt
        ORGMAT = s.ORGMAT + r.RORGMAT * delt 
        CORG = s.CORG + r.RCORG * delt 
        NORG= s.NORG + r.RNORG * delt 

        # Calculate amounts of NH4-N and NO3-N in next time step
        NH4 = s.NH4 + r.RNH4 * delt
        NO3 = s.NO3 + r.RNO3 * delt

        # Update state variables
        s.AGE = AGE
        s.ORGMAT = ORGMAT
        s.CORG = CORG
        s.NORG = NORG
        s.NH4 = NH4
        s.NO3 = NO3

        # Get external state variable
        RD_m = self.get_root_length(k)
        SM = self.get_soil_moisture_content(k)

        # Calculate the amount of N available for root uptake in the next time step
        s.NAVAIL = sinm.calculate_NAVAIL(self.soiln_profile, p.KSORP, s.NH4, s.NO3, RD_m, SM) / self.m2_to_ha
        self.check_mass_balances(day, delt)

        # Set output variables
        s.ORGMATT = np.sum(s.ORGMAT)  * (1/self.m2_to_ha)
        s.CORGT = np.sum(s.CORG)  * (1/self.m2_to_ha)
        s.NORGT = np.sum(s.NORG)  * (1/self.m2_to_ha)
        s.RMINT += np.sum(r.RNORGDIS) * (1/self.m2_to_ha)
        s.NH4T = np.sum(s.NH4) * (1/self.m2_to_ha)
        s.NO3T = np.sum(s.NO3) * (1/self.m2_to_ha)
        NH4LEACHCUM = s.NH4LEACHCUM + r.RNH4LEACHCUM * delt
        NO3LEACHCUM = s.NO3LEACHCUM + r.RNO3LEACHCUM * delt
        NDENITCUM = s.NDENITCUM + r.RNDENITCUM * delt
        NLOSSCUM =  NH4LEACHCUM + NO3LEACHCUM + NDENITCUM
        s.NH4LEACHCUM = NH4LEACHCUM
        s.NO3LEACHCUM = NO3LEACHCUM
        s.DENITCUM = NDENITCUM
        s.NLOSSCUM = NLOSSCUM


    def _on_APPLY_N(self, amount=None, application_depth = None, cnratio=None, f_orgmat=None, f_NH4N = None, f_NO3 = None, initial_age =None):
        r = self.rates
        s = self.states
        delt = 1.

        # Create model components
        samm = self.SoilInorganicNModel.SoilAmmoniumNModel()
        sni = self.SoilInorganicNModel.SoilNNitrateModel()
        sonm = self.SoilOrganicNModel()

        # Initialize ammendment rates
        RAGE_am = np.zeros((1, len(self.soiln_profile)))
        AGE0_am =  np.zeros_like(RAGE_am)
        RORGMAT_am = np.zeros_like(RAGE_am)
        RCORG_am =  np.zeros_like(RAGE_am)
        RNORG_am =  np.zeros_like(RAGE_am)
        RNH4_am = np.zeros_like(s.NH4)
        RNO3_am = np.zeros_like(s.NO3)

        # Make matrix for new amendment
        zmin = 0.
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            zmax = zmin + self.soiln_profile[il].Thickness
            AGE0_am[0,il] = initial_age * self.y_to_d
            RAGE_am[0,il] = initial_age * self.y_to_d            
            RORGMAT_am[0, il] = sonm.calculate_organic_material_application_amount(amount, application_depth, f_orgmat, layer.Thickness, zmax, zmin) * self.m2_to_ha
            RCORG_am[0, il] =  sonm.calculate_organic_carbon_application_amount(amount, application_depth, f_orgmat, layer.Thickness, zmax, zmin) * self.m2_to_ha
            RNORG_am[0, il] = sonm.calculate_organic_nitrogen_application_amount(amount, application_depth, cnratio, f_orgmat, layer.Thickness, zmax, zmin) * self.m2_to_ha
            RNH4_am[il] = samm.calculate_NH4_application_amount(amount, application_depth, f_NH4N, layer.Thickness, zmax, zmin) * self.m2_to_ha
            RNO3_am[il] = sni.calculate_NO3_application_amount(amount, application_depth, f_NO3, layer.Thickness, zmax, zmin) * self.m2_to_ha
            zmin = zmax

        # Add a new column to the state variables for organic ammendments to add the new ammendment.
        s.AGE0 = np.concatenate((s.AGE0, AGE0_am), axis = 0)
        s.ORGMAT = np.concatenate((s.ORGMAT, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.CORG = np.concatenate((s.CORG, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.NORG = np.concatenate((s.NORG, np.zeros((1, len(self.soiln_profile)))), axis = 0)
        s.AGE = np.concatenate((s.AGE, np.zeros((1, len(self.soiln_profile)))), axis = 0)

        # Store the amendment rates
        self._RAGEAM = np.concatenate((self._RAGEAM, RAGE_am), axis = 0)
        self._RORGMATAM = np.concatenate((self._RORGMATAM, RORGMAT_am), axis = 0)
        self._RCORGAM = np.concatenate((self._RCORGAM, RCORG_am), axis = 0)
        self._RNORGAM = np.concatenate(( self._RNORGAM, RNORG_am), axis = 0)
        self._RNH4AM = RNH4_am
        self._RNO3AM = RNO3_am

    def get_soil_moisture_content(self, k):
        SM = k.SM
        return SM

    def get_water_flow_rates(self, k):
        flow_m_per_d = k.Flow * self.cm_to_m
        return flow_m_per_d

    def get_N_demand(self, k):
        if "RNuptake" in k:
            N_demand_soil = k.RNuptake * self.m2_to_ha
        else:
            N_demand_soil = 0.
        return N_demand_soil

    def get_root_length(self, k):
        if "RD" in k:
            RD_m = k.RD * self.cm_to_m
        else:
            RD_m = 0.
        return RD_m

    def check_mass_balances(self, day, delt):
        s = self.states
        r = self.rates

        s.RORGMATAMTT += delt * r.RORGMATAM.sum()
        s.RORGMATDISTT += delt * r.RORGMATDIS.sum()
        s.RCORGAMTT += delt * r.RCORGAM.sum()
        s.RCORGDISTT += delt * r.RCORGDIS.sum()
        s.RNORGAMTT += delt * r.RNORGAM.sum()
        s.RNORGDISTT += delt * r.RNORGDIS.sum()

        s.RNH4MINTT+= delt * r.RNH4MIN.sum()
        s.RNH4NITRTT+= delt * r.RNH4NITR.sum()
        s.RNH4UPTT+= delt * r.RNH4UP.sum()
        s.RNH4INTT+= delt * r.RNH4IN.sum()
        s.RNH4OUTTT+= delt * r.RNH4OUT.sum()
        s.RNH4AMTT+= delt * r.RNH4AM.sum()

        s.RNO3NITRTT+=  delt * r.RNO3NITR.sum()
        s.RNO3DENITRTT+= delt * r.RNO3DENITR.sum()
        s.RNO3UPTT+= delt * r.RNO3UP.sum()
        s.RNO3INTT+= delt * r.RNO3IN.sum()
        s.RNO3OUTTT+= delt * r.RNO3OUT.sum()
        s.RNO3AMTT+= delt * r.RNO3AM.sum()

        ORGMATBAL = self._ORGMATI.sum() - s.ORGMAT.sum() + s.RORGMATAMTT - s.RORGMATDISTT
        if(abs(ORGMATBAL) > 0.0001):
            msg = "Organic matter mass balance is not closing on %s with checksum: %f" % (day, ORGMATBAL)
            raise exc.SoilOrganicMatterBalanceError(msg)

        CORGBAL = self._CORGI.sum() - s.CORG.sum() + s.RCORGAMTT - s.RCORGDISTT
        if(abs(CORGBAL) > 0.0001):
            msg = "Organic carbon mass balance is not closing on %s with checksum: %f" % (day, CORGBAL)
            raise exc.SoilOrganicCarbonBalanceError(msg)

        NORGBAL = self._NORGI.sum() - s.NORG.sum() + s.RNORGAMTT - s.RNORGDISTT
        if(abs(NORGBAL) > 0.0001):
            msg = "Organic carbon mass balance is not closing on %s with checksum: %f" % (day, NORGBAL)
            raise exc.SoilOrganicNitrogenBalanceError(msg)

        NH4BAL = self._NH4I.sum() - s.NH4.sum() + s.RNH4AMTT + s.RNH4INTT + s.RNH4MINTT - s.RNH4NITRTT - s.RNH4OUTTT - s.RNH4UPTT
        if(abs(NH4BAL) > 0.0001):
            msg = "NH4-N mass balance is not closing on %s with checksum: %f" % (day, NH4BAL)
            raise exc.SoilAmmoniumBalanceError(msg)

        NO3BAL = self._NO3I.sum() - s.NO3.sum() + s.RNO3AMTT + s.RNO3NITRTT + s.RNO3INTT - s.RNO3DENITRTT - s.RNO3OUTTT - s.RNO3UPTT
        if(abs(NO3BAL) > 0.0001):
            msg = "NO3-N mass balance is not closing on %s with checksum: %f" % (day, NO3BAL)
            raise exc.SoilAmmoniumBalanceError(msg)

    class SoilInorganicNModel():
        def calculate_flow_rates(self, soiln_profile, flow_m_per_d, KSORP, NH4, NO3, SM):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4IN, RNH4OUT = samm.calculate_NH4_flow_rates(soiln_profile, flow_m_per_d, KSORP, NH4, SM)
            RNO3IN, RNO3OUT = sni.calculate_NO3_flow_rates(soiln_profile, flow_m_per_d, NO3, SM)
            return RNH4IN, RNH4OUT, RNO3IN, RNO3OUT

        def calculate_reaction_rates(self, soiln_profile, KDENIT_REF, KNIT_REF, KSORP, MRCDIS, NH4, NO3, RCORGDIS, RNORGDIS, SM, T, WFPS_CRIT):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            RNH4MIN, RNH4NITR = samm.calculate_NH4_reaction_rates(soiln_profile, KNIT_REF, KSORP, NH4, RNORGDIS, SM, T)
            RNO3NITR, RNO3DENITR = sni.calculate_NO3_reaction_rates(soiln_profile, KDENIT_REF, MRCDIS, NO3, RCORGDIS, RNH4NITR, SM, T, WFPS_CRIT)  
            return RNH4MIN, RNH4NITR, RNO3NITR, RNO3DENITR

        def calculate_NAVAIL(self, soiln_profile, KSORP, NH4, NO3, RD_m, SM):
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            zmin = 0.
            NAVAIL = 0.
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + layer.Thickness_m
                NH4_avail_layer = samm.calculate_available_NH4(KSORP, NH4[il], RD_m, layer.RHOD_kg_per_m3, SM[il], zmax, zmin)
                NO3_avail_layer = sni.calculate_available_NO3(NO3[il], RD_m, zmax, zmin)     
                NAVAIL += (NH4_avail_layer + NO3_avail_layer)
                zmin = zmax
            return NAVAIL

        def calculate_N_uptake_rates(self, soiln_profile, delt, KSORP, N_demand_soil, NH4, NO3, RD_m, SM):
            RNH4UP = np.zeros_like(NH4)
            RNO3UP = np.zeros_like(NO3)
            samm = self.SoilAmmoniumNModel()
            sni = self.SoilNNitrateModel()
            zmin = 0.
            for il, layer in enumerate(soiln_profile):
                zmax = zmin + layer.Thickness_m
                RNH4UP[il] = samm.calculate_NH4_plant_uptake_rate(KSORP, N_demand_soil, NH4[il], RD_m, layer.RHOD_kg_per_m3, SM[il], zmax, zmin)
                N_demand_soil-= RNH4UP[il] * delt
                RNO3UP[il] = sni.calculate_NO3_plant_uptake_rate(N_demand_soil, NO3[il], RD_m, zmax, zmin)
                N_demand_soil-= RNO3UP[il] * delt
                zmin = zmax
            return RNH4UP, RNO3UP  

        class SoilAmmoniumNModel():
            def calculate_NH4_reaction_rates(self, soiln_profile, KNITREF, KSORP, NH4, RNORGDIS, SM, T):
                RNH4MIN = np.zeros_like(NH4)
                RNH4NITR = np.zeros_like(NH4)
                for il, layer in enumerate(soiln_profile):
                    RNMIN_kg_per_m2 = RNORGDIS[:,il].sum()
                    RNH4MIN[il] = self.calculate_mineralization_rate(RNMIN_kg_per_m2)
                    RNH4NITR[il] = self.calculate_nitrification_rate(KNITREF, KSORP, layer.Thickness_m, NH4[il], layer.RHOD_kg_per_m3, SM[il], layer.SM0, T)
                return RNH4MIN, RNH4NITR

            def calculate_NH4_flow_rates(self, soiln_profile, flow_m_per_d, KSORP, NH4, SM):
                RNH4IN = np.zeros_like(NH4)
                RNH4OUT = np.zeros_like(NH4)
                for il in range(0, len(NH4)):
                    layer = soiln_profile[il]
                    RNH4OUT[il] = self.calculate_NH4_outflow_rate(flow_m_per_d[il+1], il, KSORP, NH4[il], layer.Thickness_m, layer.RHOD_kg_per_m3, SM[il])
                    RNH4IN[il] = self.calculate_NH4_inflow_rate(il, RNH4OUT[il-1])
                return RNH4IN, RNH4OUT

            def calculate_NH4_application_amount(self, amount, application_depth, f_NH4N, layer_thickness, zmax, zmin):
                if(application_depth > zmax):
                    NH4_am = (layer_thickness / application_depth) * f_NH4N *  amount
                elif(application_depth >= zmin and application_depth <= zmax):
                    NH4_am = ((application_depth - zmin) / application_depth) * f_NH4N  * amount
                else:
                    NH4_am = 0.
                return NH4_am

            def calculate_NH4_concentration(self, KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM):
                cNH4 = (SM / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4 / (layer_thickness * SM)
                return cNH4

            def calculate_available_NH4(self, KSORP, NH4, RD, RHOD_kg_per_m3, SM, zmax, zmin):
                dz = zmax - zmin
                if(RD <= zmin):
                    NH4_avail = 0.
                elif(RD > zmax):
                    NH4_avail = (SM / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4
                else:
                    NH4_avail = ((RD - zmin)/ dz) * (SM / ( KSORP * RHOD_kg_per_m3 + SM)) * NH4
                return NH4_avail

            def calculate_mineralization_rate(self, rNMINs_layer):
                RNH4MIN = (rNMINs_layer).sum()
                return RNH4MIN

            def calculate_nitrification_rate(self, KNIT_REF, KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM, SM0, T):
                cNH4 = self.calculate_NH4_concentration(KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM)
                fWNIT = self.calculate_soil_moisture_response_nitrification_rate_constant(SM, SM0)
                fT = self.calculate_temperature_response_nitrification_rate_constant(T)
                RNH4NIT = fWNIT * fT *  KNIT_REF * SM * cNH4 * layer_thickness
                return RNH4NIT

            def calculate_NH4_plant_uptake_rate(self, KSORP, N_demand_soil, NH4, RD_m, RHOD_kg_per_m3, SM, zmax, zmin):
                NH4_av = self.calculate_available_NH4(KSORP, NH4, RD_m, RHOD_kg_per_m3, SM, zmax, zmin)
                RNH4UP = min(N_demand_soil, NH4_av)
                return RNH4UP

            def calculate_NH4_inflow_rate(self, il, RNH4OUT_above):
                if(il == 0):
                    RNH4IN = 0.
                else:               
                    RNH4IN = RNH4OUT_above
                return RNH4IN

            def calculate_NH4_outflow_rate(self, flow_m_per_d, il, KSORP, NH4, layer_thickness, RHOD_kg_per_m3, SM):
                cNH4 = self.calculate_NH4_concentration(KSORP, layer_thickness, NH4, RHOD_kg_per_m3, SM)
                if(il == 0):
                    RNH4OUT = cNH4 * max(0,  flow_m_per_d)
                else:               
                    RNH4OUT = cNH4 * flow_m_per_d
                return RNH4OUT

            def calculate_soil_moisture_response_nitrification_rate_constant(self, SM, SM0):
                WFPS = SM / SM0
                fWNIT = 0.9 / (1. + np.exp(-15 *(WFPS - 0.45))) + 0.1 - 1/(1+np.exp(-50. * (WFPS - 0.95)))
                return fWNIT

            def calculate_temperature_response_nitrification_rate_constant(self, T):
                fT = 1/(1+np.exp(-0.26*(T-17.)))-1/(1+np.exp(-0.77*(T-41.9)))
                return fT

        class SoilNNitrateModel():
            def calculate_NO3_flow_rates(self, soiln_profile, flow_m_per_d, NO3, SM):
                RNO3OUT = np.zeros_like(NO3)
                RNO3IN = np.zeros_like(NO3)
                for il in range(0, len(NO3)):
                    layer = soiln_profile[il]
                    RNO3OUT[il] = self.calculate_NO3_outflow_rate(il, flow_m_per_d[il+1], layer.Thickness_m, NO3[il], SM[il])
                    RNO3IN[il] = self.calculate_NO3_inflow_rate(il, RNO3OUT[il-1])
                return RNO3IN, RNO3OUT

            def calculate_NO3_reaction_rates(self, soiln_profile, KDENIT_REF, MRCDIS, NO3, RCORGDIS, RNH4NITR, SM, T, WFPS_CRIT):
                RNO3NITR = np.zeros_like(NO3)     
                RNO3DENITR = np.zeros_like(NO3)     
                for il, layer in enumerate(soiln_profile):           
                    RNO3NITR[il] = RNH4NITR[il]
                    RNO3DENITR[il] = self.calculate_denitrification_rate(layer.Thickness_m, NO3[il], KDENIT_REF, MRCDIS, RCORGDIS[:,il].sum(), SM[il], layer.SM0, T, WFPS_CRIT)
                return RNO3NITR, RNO3DENITR      

            def calculate_NO3_application_amount(self, amount, application_depth, f_NO3N, layer_thickness, zmax, zmin):
                if(application_depth > zmax):
                    NO3_am = (layer_thickness / application_depth) * f_NO3N *  amount
                elif(application_depth >= zmin and application_depth <= zmax):
                    NO3_am = ((application_depth - zmin) / application_depth) * f_NO3N  * amount
                else:
                    NO3_am = 0.
                return NO3_am

            def calculate_NO3_concentration(self, layer_thickness, NO3, SM):
                cNO3 = NO3 / (layer_thickness * SM)
                return cNO3

            def calculate_available_NO3(self, NO3, RD, zmax, zmin):
                dz = zmax - zmin
                if(RD <= zmin):
                    NO3_avail_layer = 0.
                elif(RD > zmax):
                    NO3_avail_layer = NO3            
                else:
                    NO3_avail_layer = ((RD - zmin)/ dz) * NO3 
                return NO3_avail_layer

            def calculate_denitrification_rate(self, layer_thickness, NO3, KDENIT_REF, MRCDIS, RCORGT_kg_per_m2, SM, SM0, T, WFPS_CRIT):
                cNO3 = self.calculate_NO3_concentration(layer_thickness, NO3, SM)
                fR = self.calculate_soil_respiration_response_denitrifiation_rate_constant(RCORGT_kg_per_m2, MRCDIS)
                fW = self.calculate_soil_moisture_response_denitrification_rate_constant(SM, SM0, WFPS_CRIT)
                fT = self.calculate_temperature_response_denitrification_rate_constant(T)
                RNO3DENIT = fW * fT * fR * KDENIT_REF * SM * cNO3 * layer_thickness
                return RNO3DENIT

            def calculate_NO3_plant_uptake_rate(self, N_demand_soil, NO3, RD_m, zmax, zmin):
                NO3_av = self.calculate_available_NO3(NO3, RD_m, zmax, zmin)
                NO3UPT_kg_per_m2 = min(N_demand_soil, NO3_av)
                N_demand_soil -= NO3UPT_kg_per_m2
                RNO3UP = NO3UPT_kg_per_m2
                return RNO3UP

            def calculate_NO3_inflow_rate(self, il, RNO3OUT_above):
                if(il == 0):
                    RNO3IN = 0.
                else:               
                    RNO3IN = RNO3OUT_above
                return RNO3IN

            def calculate_NO3_outflow_rate(self, il, flow_m_per_d, layer_thickness, NO3, SM):
                cNO3 = self.calculate_NO3_concentration(layer_thickness, NO3, SM)
                if(il == 0):
                    RNO3OUT = cNO3 * max(0,  flow_m_per_d)
                else:               
                    RNO3OUT = cNO3 * flow_m_per_d
                return RNO3OUT

            def calculate_soil_moisture_response_denitrification_rate_constant(self, SM, SM0, WFPS_CRIT):
                WFPS = SM / SM0
                if(WFPS < WFPS_CRIT):
                    fW = 0.
                else:
                    fW = np.power((WFPS - WFPS_CRIT)/(1 - WFPS_CRIT),2)
                return fW

            def calculate_soil_respiration_response_denitrifiation_rate_constant(self, RCORGT, MRCDIS):
                fR = RCORGT / (MRCDIS + RCORGT)
                return fR

            def calculate_temperature_response_denitrification_rate_constant(self, T):
                fT = 1/(1+np.exp(-0.26*(T-17)))-1/(1+np.exp(-0.77*(T-41.9)))
                return fT

    class SoilOrganicNModel():
        def calculate_apparent_age_increase_rate(self, AGE, delt, T):
            RAGEAG = np.zeros_like(AGE)
            janssen = self.Janssen()
            for am in range(0, AGE.shape[0]):
                for il in range(0, AGE.shape[1]):
                    RAGEAG[am,il] = janssen.calculate_increase_apparent_age_rate(delt, T)
            return RAGEAG

        def calculate_dissimilation_rates(self, AGE, AGE0, CNRatioBio, FASDIS, NORG, ORGMAT, T):
            RORGMATDIS = np.zeros_like(AGE)
            RCORGDIS = np.zeros_like(AGE)
            RNORGDIS = np.zeros_like(AGE)
            janssen = self.Janssen()
            minip_c = self.MINIP_C()
            minip_n = self.MINIP_N()

            for am in range(0, AGE.shape[0]):
                for il in range(0, AGE.shape[1]):
                    #r.RAGE[am,il] = janssen.calculate_increase_apparent_age_rate(delt, T)
                    if(ORGMAT[am, il] > 0):
                        RORGMATDIS[am,il] = janssen.calculate_dissimilation_rate_OM(ORGMAT[am,il], AGE0[am,il], AGE[am,il])
                        RCORGDIS[am,il] = minip_c.calculate_dissimilation_rate_C(janssen, ORGMAT[am,il], AGE0[am,il], AGE[am,il], T)
                        RNORGDIS[am,il] = minip_n.calculate_dissimilation_rate_N(janssen, minip_c, ORGMAT[am,il], NORG[am,il], AGE0[am,il], FASDIS, CNRatioBio, AGE[am,il], T)
                    else:
                        RORGMATDIS[am,il] = 0.
                        RCORGDIS[am,il] = 0.                    
                        RNORGDIS[am,il] = 0.
            return RORGMATDIS, RCORGDIS, RNORGDIS

        def calculate_organic_carbon_application_amount(self, amount, application_depth, f_orgmat, layer_thickness, zmax, zmin):
            minip_C = self.MINIP_C()
            ORGMAT_am = self.calculate_organic_material_application_amount(amount, application_depth, f_orgmat, layer_thickness, zmax, zmin)
            CORG_am = minip_C.calculate_organic_C(ORGMAT_am)
            return CORG_am

        def calculate_organic_nitrogen_application_amount(self, amount, application_depth, cnratio, f_orgmat, layer_thickness, zmax, zmin):
            minip_C = self.MINIP_C()
            ORGMAT_am = self.calculate_organic_material_application_amount(amount, application_depth, f_orgmat, layer_thickness, zmax, zmin)
            CORG_am = minip_C.calculate_organic_C(ORGMAT_am)
            if(cnratio == 0):
                NORG_am = 0.
            else:
                NORG_am = CORG_am / cnratio
            return NORG_am

        def calculate_organic_material_application_amount(self, amount, application_depth, f_orgmat, layer_thickness, zmax, zmin):
            if(application_depth > zmax):
                ORGMAT_am = (layer_thickness / application_depth) * f_orgmat * amount
            elif(application_depth >= zmin and application_depth <= zmax):
                ORGMAT_am = ((application_depth - zmin) / application_depth) * f_orgmat * amount
            else:
                ORGMAT_am = 0
            return ORGMAT_am

        class Janssen():
            m = 1.6
            b = 2.82
            y_to_d = 365.25

            def calculate_increase_apparent_age_rate(self, dt, T):
                f_T = self.calculate_temperature_response_dissimilation_rate(T)
                dA = f_T * dt
                return dA

            def calculate_organic_matter_amount_analytical(self, a, OM0, t):
                m = self.m
                b = self.b        
                OM = OM0 * np.exp((b/(m-1)) * (pow((a + t)/self.y_to_d,(1-m)) - pow((a/self.y_to_d),1-m)))
                return OM

            def calculate_relative_dissimilation_rate_OM(self, a, t):
                m = self.m
                b = self.b   
                k = b * pow((a + t)/self.y_to_d, -m) / self.y_to_d
                return k

            def calculate_relative_dissimilation_rate_OM_T(self, a, t, T):
                m = self.m
                b = self.b   
                f_T = self.calculate_temperature_response_dissimilation_rate(T)
                k = b  * f_T * pow((a + t * f_T)/self.y_to_d, -m) / self.y_to_d
                return k

            def calculate_dissimilation_rate_OM(self, OM, a, t): 
                k = self.calculate_relative_dissimilation_rate_OM(a, t)
                rate = k * OM
                return rate

            def calculate_dissimilation_rate_OM_T(self, OM, a, t, T):
                k = self.calculate_relative_dissimilation_rate_OM_T(a, t, T)
                rate = k * OM
                return rate

            def calculate_temperature_response_dissimilation_rate(self, T):
                f_T = pow(2, (T-9)/9)
                return f_T

        class MINIP_C(object):
            OM_to_C = 0.58
            y_to_d = 365.

            def calculate_assimilation_rate(self, janssen, OM, a, f_ass_dis, t, T):
                r_disc = self.calculate_dissimilation_rate_C(janssen, OM, a, t, T)
                r_ass = r_disc * f_ass_dis
                return r_ass
    
            def calculate_dissimilation_rate_C(self, janssen, OM, a, t, T):
                k = janssen.calculate_relative_dissimilation_rate_OM_T(a, t, T)
                Corg = self.calculate_organic_C(OM)
                rate = k * Corg
                return rate

            def calculate_total_conversion_rate_C(self, janssen, OM, a, f_ass_dis, t, T):
                r_dis_C = self.calculate_dissimilation_rate_C(janssen, OM, a, t, T)
                r_ass_C = self.calculate_assimilation_rate(janssen, OM, a, f_ass_dis, t, T)
                r_conv_C = r_dis_C + r_ass_C
                return r_conv_C

            def calculate_organic_C(self, OM):
                Corg = OM * self.OM_to_C
                return Corg

        class MINIP_N(object):    

            def calculate_total_conversion_rate_N(self, janssen, minip_c, OM, Norg, a, f_ass_dis, t, T):
                r_conv_C = minip_c.calculate_total_conversion_rate_C(janssen, OM, a, f_ass_dis, t, T)
                C = minip_c.calculate_organic_C(OM)
                r_conv_N = r_conv_C * (Norg/C)
                return r_conv_N

            def calculate_assimilation_rate_N(self, janssen, minip_c, OM, a, f_ass_dis, f_C_N_microbial, t, T):
                r_ass_C = minip_c.calculate_assimilation_rate(janssen, OM, a, f_ass_dis, t, T)
                r_ass_N = r_ass_C/f_C_N_microbial
                return r_ass_N

            def calculate_dissimilation_rate_N(self, janssen, minip_c, OM, Norg, a, f_ass_dis, f_C_N_microbial, t, T):
                r_ass_N = self.calculate_assimilation_rate_N(janssen, minip_c, OM, a, f_ass_dis, f_C_N_microbial, t, T)
                r_conv_N = self.calculate_total_conversion_rate_N(janssen, minip_c, OM, Norg, a, f_ass_dis, t, T)
                r_diss_N = r_conv_N - r_ass_N
                return r_diss_N
