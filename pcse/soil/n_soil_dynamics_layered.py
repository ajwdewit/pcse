# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf
import numpy as np
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
        
        RNH4 = Instance(np.ndarray)
        RNH4MIN = Instance(np.ndarray)
        RNH4NITR = Instance(np.ndarray)
        RNH4UP = Instance(np.ndarray)
        RNH4IN = Instance(np.ndarray)
        RNH4OUT = Instance(np.ndarray)

        RNO3 = Instance(np.ndarray)
        RNO3NITR = Instance(np.ndarray)
        RNO3DENITR = Instance(np.ndarray)
        RNO3UP = Instance(np.ndarray)
        RNO3IN = Instance(np.ndarray)
        RNO3OUT = Instance(np.ndarray)

        RNH4LEACHCUM = Float()
        RNO3LEACHCUM = Float()
        RNDENITCUM = Float()
        RNLOSS = Float()

    class Parameters(ParamTemplate):
        A0SOM = Float()            # Initial age of humus (y)
        CNRatioBio = Float()       # C:N Ratio of microbial biomass (kg C kg-1 N)
        FASDIS = Float()           # Fraction of assimilation to dissimilation (kg kg-1)
        KDENIT_REF = Float()       # Reference first order denitrification rate constant (d-1)
        KNIT_REF = Float()         # Reference first order nitrification rate constant (d-1)
        KSORP = Float()            # Sorption coefficient ammonium (m3 kg-1)
        MRCDIS = Float()           # Michaelis Menten constant for response factor denitrification to soil respiration
        WFPS_CRIT = Float()        # Critical water filled pore space fraction (m3 water m-3 pore) for denitrification

    def initialize(self, day, kiosk, parvalues):
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)
        self.soiln_profile = SoilNProfile(parvalues)

        # Initialize states
        parvalues._soildata["soil_profile"] = self.soiln_profile
        NH4 = np.zeros(len(self.soiln_profile))
        NO3 = np.zeros(len(self.soiln_profile))
        AGE =  np.zeros((1, len(self.soiln_profile)))
        AGE0 = np.zeros_like(AGE)
        ORGMAT = np.zeros_like(AGE)
        CORG =  np.zeros_like(AGE)
        NORG =  np.zeros_like(AGE)
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            NH4[il] = layer.NH4I
            NO3[il] = layer.NO3I
            AGE0[0,il] = self.params.A0SOM * self.y_to_d
            AGE[0,il] = self.params.A0SOM * self.y_to_d
            ORGMAT[0,il] = layer.RHOD * layer.FSOMI * layer.Thickness *  self.g_to_kg / self.cm2_to_ha
            CORG[0,il] = minip_C.calculate_organic_C(ORGMAT[0,il])
            NORG[0,il] = CORG[0, il] / layer.CNRatioSOMI

        CORGT = np.sum(CORG)
        NORGT = np.sum(NORG)
        ORGMATT = np.sum(ORGMAT)
        RMINT = 0.
        NH4T = np.sum(NH4)
        NO3T = np.sum(NO3)
        NAVAIL = 0.
        NH4LEACHCUM = 0.
        NO3LEACHCUM = 0.
        NDENITCUM = 0.
        NLOSSCUM = 0.

        samm = self.SoilAmmoniumNModel()
        sni = self.SoilNNitrateModel()

        zmin = 0.
        for il, layer in enumerate(self.soiln_profile):
            dz = layer.Thickness
            zmax = zmin + dz
            RHOD = layer.RHOD
            RHOD_kg_per_m2 = RHOD * self.g_to_kg / self.cm3_to_m3
            NH4_avail_layer = samm.calculate_available_NH4(self.params.KSORP, NH4[il], 0., RHOD_kg_per_m2, self.kiosk.SM[il], zmax, zmin)
            NO3_avail_layer = sni.calculate_available_NO3(NO3[il], 0., zmax, zmin)     
            NAVAIL += NH4_avail_layer + NO3_avail_layer
            zmin = zmax

        states = {"NAVAIL": NAVAIL, "NO3": NO3, "NH4": NH4, "AGE": AGE, "AGE0": AGE0, "ORGMAT": ORGMAT, "CORG": CORG, "NORG": NORG,  
                  "ORGMATT": ORGMATT,  "CORGT": CORGT, "NORGT": NORGT, "RMINT": RMINT, "NH4T": NH4T, "NO3T": NO3T, 
                  "NH4LEACHCUM": NH4LEACHCUM, "NO3LEACHCUM": NO3LEACHCUM, "NDENITCUM": NDENITCUM, "NLOSSCUM": NLOSSCUM}
        #self.states = self.StateVariables(kiosk, publish=["NAVAIL", "ORGMAT", "CORG", "NORG", "ORGMATT", "CORGT", "NORGT"], **states)
        self.states = self.StateVariables(kiosk, publish=["NAVAIL", "ORGMATT", "CORGT", "NORGT"], **states)
        self.rates = self.RateVariables(kiosk)

        self._connect_signal(self._on_APPLY_N, signals.apply_n)

    @prepare_rates
    def calc_rates(self, day, drv):
        k = self.kiosk
        r = self.rates
        s = self.states
        p = self.params

        delt = 1.0
        T = drv.TEMP

        janssen = self.SoilOrganicNModel.Janssen()
        minip_c = self.SoilOrganicNModel.MINIP_C()
        minip_n = self.SoilOrganicNModel.MINIP_N()

        # initialize rates
        r.RAGE = np.zeros((1, len(self.soiln_profile)))
        r.RORGMAT = np.zeros_like(r.RAGE)
        r.RCORG = np.zeros_like(r.RAGE)
        r.RNORG = np.zeros_like(r.RAGE)

        for am in range(0, r.RAGE.shape[0]):
            for il in range(0, r.RAGE.shape[1]):
                r.RAGE[am,il] = janssen.calculate_increase_apparent_age_rate(delt, T)
                if(s.ORGMAT[am, il] > 0):
                    r.RORGMAT[am,il] = janssen.calculate_dissimilation_rate_OM(s.ORGMAT[am,il], s.AGE0[am,il], s.AGE[am,il])
                    r.RCORG[am,il] = minip_c.calculate_dissimilation_rate_C(janssen, s.ORGMAT[am,il], s.AGE0[am,il], s.AGE[am,il], T)
                    r.RNORG[am,il] = minip_n.calculate_dissimilation_rate_N(janssen, minip_c, s.ORGMAT[am,il], s.NORG[am,il], s.AGE0[am,il], p.FASDIS, p.CNRatioBio, s.AGE[am,il], T)
                else:
                    r.RORGMAT[am,il] = 0.
                    r.RCORG[am,il] = 0.
                    r.RNORG[am,il] = 0.

        # initialize rates for ammonium
        r.RNH4 = np.zeros(len(self.soiln_profile))
        r.RNH4MIN = np.zeros_like(len(r.RNH4 )
        r.RNH4NITR = np.zeros_like(len(r.RNH4 )
        r.RNH4UP = np.zeros_like(len(r.RNH4 )
        r.RNH4IN = np.zeros_like(len(r.RNH4 )
        r.RNH4OUT = np.zeros_like(len(r.RNH4 )

        # Initialize rates for nitrate
        r.RNO3  = np.zeros_like(len(r.RNH4 )
        r.RNO3NITR = np.zeros_like(len(r.RNH4 )
        r.RNO3DENITR = np.zeros_like(len(r.RNH4 )
        r.RNO3UP = np.zeros_like(len(r.RNH4 )
        r.RNO3IN = np.zeros_like(len(r.RNH4 )
        r.RNO3OUT = np.zeros_like(len(r.RNH4 )
        r.RDENITCUM = 0.

        # Calculate N uptake
        samm = self.SoilAmmoniumNModel()
        sni = self.SoilNNitrateModel()

        if "RNuptake" in k:
            N_demand_soil = k.RNuptake
        else:
            N_demand_soil = 0.
        if "RD" in k:
            RD = k.RD
        else:
            RD = 0.

        NH4UPT_kg_per_ha = np.zeros(len(self.soiln_profile))
        NO3UPT_kg_per_ha = np.zeros(len(self.soiln_profile))

        zmin = 0.
        for il, layer in enumerate(self.soiln_profile):
            dz = layer.Thickness 
            zmax = zmin + dz
            RHOD = layer.RHOD
            RHOD_kg_per_m2 = RHOD * self.g_to_kg / self.cm3_to_m3
            NH4_av = samm.calculate_available_NH4(p.KSORP, s.NH4[il], RD, RHOD_kg_per_m2, k.SM[il], zmax, zmin)
            NH4UPT_kg_per_ha[il] = min(N_demand_soil, NH4_av)
            r.RNH4UP[il] =  self.ha_to_m2 *  NH4UPT_kg_per_ha[il] / (dz * self.cm_to_m)
            N_demand_soil -= NH4UPT_kg_per_ha[il]
            NO3_av = sni.calculate_available_NO3(s.NO3[il], RD, zmax, zmin)
            NO3UPT_kg_per_ha[il] = min(N_demand_soil, NO3_av)
            N_demand_soil -= NO3UPT_kg_per_ha[il] 
            r.RNO3UP[il] = self.ha_to_m2 *  NO3UPT_kg_per_ha[il] / (dz * self.cm_to_m)

        # Calculate reactions ammonium and nitrate
        samm = self.SoilAmmoniumNModel()
        sni = self.SoilNNitrateModel()
        NO3PRE = s.NO3 + r.RNO3UP * delt
        NH4PRE = s.NH4 + r.RNH4UP * delt

        for il, layer in enumerate(self.soiln_profile):
            zmax = zmin + dz
            RHOD = layer.RHOD
            RHOD_kg_per_m2 = RHOD * self.g_to_kg / self.cm3_to_m3
            dz = self.soiln_profile[il].Thickness * self.cm_to_m
            SM0 = self.soiln_profile[il].SM0
            RNMIN_kg_per_m3 = r.RNORG[:,il] * self.m2_to_ha

            cNH4 = (k.SM[il] / ( p.KSORP * RHOD_kg_per_m2 + k.SM[il])) * NH4PRE[il] * self.m2_to_ha  / (dz * k.SM[il])
            r.RNH4MIN[il] = samm.calculate_mineralization_rate(dz, RNMIN_kg_per_m3)
            r.RNH4NITR[il] = samm.calculate_nitrification_rate(cNH4, p.KNIT_REF, k.SM[il],SM0, T)
            r.RNH4[il] = (1/self.m2_to_ha) * dz * (r.RNH4MIN[il] - r.RNH4NITR[il] - r.RNH4UP[il])

            cNO3 = NO3PRE[il] * self.m2_to_ha / (dz * k.SM[il])
            RCORGT_kg_per_m2 = - r.RCORG.sum() * self.ha_to_m2
            SM0 = self.soiln_profile[il].SM0
            r.RNO3NITR[il] = r.RNH4NITR[il]
            r.RNO3DENITR[il] = sni.calculate_denitrification_rate(cNO3, p.KDENIT_REF, p.MRCDIS, RCORGT_kg_per_m2, k.SM[il], SM0, T, p.WFPS_CRIT)
            r.RDENITCUM += (1/self.m2_to_ha) * dz * r.RNO3DENITR[il]
            r.RNO3[il] =  (1/self.m2_to_ha) * dz *  (r.RNO3NITR[il]  - r.RNO3DENITR[il] - r.RNO3UP[il])

        NH4PRE2 = s.NH4  + r.RNH4 * delt
        NO3PRE2 = s.NO3 + r.RNO3 * delt

        # Calculate flow rates
        flow_m_per_d = k.Flow * self.cm_to_m

        for il in range(0, len(NO3PRE)):
            layer = self.soiln_profile[il] 
            zmax = zmin + dz
            RHOD = layer.RHOD
            RHOD_kg_per_m2 = RHOD * self.g_to_kg / self.cm3_to_m3
            dz = self.soiln_profile[il].Thickness * self.cm_to_m
            cNH4 = (k.SM[il] / ( p.KSORP * RHOD_kg_per_m2 + k.SM[il])) * NH4PRE2[il] * self.m2_to_ha  / (dz * k.SM[il])
            cNO3 = NO3PRE2[il] * self.m2_to_ha / (dz * k.SM[il])

            if(il == 0):
                r.RNH4IN[il] = 0.
                r.RNO3IN[il] = 0.
                r.RNH4OUT[il] = cNH4 * max(0,  k.Flow[1]) * self.cm_to_m / dz
                r.RNO3OUT[il] = cNO3 * max(0,  k.Flow[1]) * self.cm_to_m / dz
            else:               
                r.RNH4IN[il] = r.RNH4OUT[il-1]
                r.RNO3IN[il] = r.RNO3OUT[il-1]
                r.RNH4OUT[il] = cNH4 * k.Flow[il+1] * self.cm_to_m / dz
                r.RNO3OUT[il] = cNO3 * k.Flow[il+1] * self.cm_to_m / dz

            r.RNO3[il] =  (1/self.m2_to_ha) * dz * (r.RNO3NITR[il]  - r.RNO3DENITR[il] - r.RNO3UP[il] + r.RNO3IN[il] - r.RNO3OUT[il])
            r.RNH4[il] =  (1/self.m2_to_ha) * dz * (r.RNH4MIN[il] - r.RNH4NITR[il] - r.RNH4UP[il] + r.RNH4IN[il] - r.RNH4OUT[il])

        r.RNH4LEACHCUM =  self.cm_to_m * self.soiln_profile[-1].Thickness * (1/self.m2_to_ha) * r.RNH4OUT[-1]
        r.RNO3LEACHCUM =  self.cm_to_m * self.soiln_profile[-1].Thickness * (1/self.m2_to_ha) * r.RNO3OUT[-1]

 
    @prepare_states
    def integrate(self, day, delt=1.0):
        k = self.kiosk
        r = self.rates
        s = self.states
        p = self.params

        s.NAVAIL = s.NAVAIL

        AGE =  np.zeros((1, len(self.soiln_profile)))
        ORGMAT = np.zeros_like(AGE)
        CORG = np.zeros_like(AGE)
        NORG = np.zeros_like(AGE)
        NO3 = np.zeros(len(self.soiln_profile))
        NH4 = np.zeros_like(NO3)

        for am in range(0, r.RAGE.shape[0]):
            for il in range(0, r.RAGE.shape[1]):
                AGE[am, il] = s.AGE[am,il] + r.RAGE[am, il] * delt
                ORGMAT[am, il] = s.ORGMAT[am,il] + r.RORGMAT[am, il] * delt
                CORG[am, il] = s.CORG[am,il] + r.RCORG[am, il] * delt
                NORG[am, il] = s.NORG[am, il] + r.RNORG[am, il] * delt

        for il in range(0, len(s.NH4)):
            NH4[il] = s.NH4[il] + r.RNH4[il] * delt
            NO3[il] = s.NO3[il] + r.RNO3[il] * delt

        NH4LEACHCUM = s.NH4LEACHCUM + r.RNH4LEACHCUM * delt
        NO3LEACHCUM = s.NO3LEACHCUM + r.RNO3LEACHCUM * delt
        NDENITCUM = s.NDENITCUM + r.RNDENITCUM * delt
        NLOSSCUM =  NH4LEACHCUM + NO3LEACHCUM + NDENITCUM

        s.AGE = AGE
        s.ORGMAT = ORGMAT
        s.CORG = CORG
        s.NORG = NORG
        s.NH4 = NH4
        s.NO3 = NO3

        s.ORGMATT = np.sum(s.ORGMAT)
        s.CORGT = np.sum(s.CORG)
        s.NORGT = np.sum(s.NORG)
        s.RMINT += np.sum(r.RNORG) * delt
        s.NH4T = np.sum(s.NH4)
        s.NO3T = np.sum(s.NO3)
        s.NH4LEACHCUM = NH4LEACHCUM
        s.NO3LEACHCUM = NO3LEACHCUM
        s.DENITCUM = NDENITCUM
        s.NLOSSCUM = NLOSSCUM

        NAVAIL = 0.

        samm = self.SoilAmmoniumNModel()
        sni = self.SoilNNitrateModel()

        if "RD" in self.kiosk:
            RD = self.kiosk.RD
        else:
            RD = 0.
        zmin = 0.
        for il, layer in enumerate(self.soiln_profile):
            dz = layer.Thickness
            zmax = zmin + dz
            RHOD = layer.RHOD
            RHOD_kg_per_m2 = RHOD * self.g_to_kg / self.cm3_to_m3
            NH4_avail_layer = samm.calculate_available_NH4(self.params.KSORP, NH4[il], RD, RHOD_kg_per_m2, self.kiosk.SM[il], zmax, zmin)
            NO3_avail_layer = sni.calculate_available_NO3(NO3[il], RD, zmax, zmin)     
            NAVAIL += NH4_avail_layer + NO3_avail_layer
            zmin = zmax

        s.NAVAIL = NAVAIL

    def _on_APPLY_N(self, amount=None, application_depth = None, cnratio=None, f_orgmat=None, f_NH4N = None, f_NO3 = None, initial_age =None):
        r = self.rates
        s = self.states

        AGE_am = np.zeros((1, len(self.soiln_profile)))
        AGE0_am =  np.zeros_like(AGE_am)
        ORGMAT_am = np.zeros_like(AGE_am)
        CORG_am =  np.zeros_like(AGE_am)
        NORG_am =  np.zeros_like(AGE_am)

        # Make matrix for new amendment
        zmin = 0.
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            AGE0_am[0,il] = initial_age * self.y_to_d
            AGE_am[0,il] = initial_age * self.y_to_d            
            zmax = zmin + self.soiln_profile[il].Thickness
            if(application_depth > zmax):
                ORGMAT_am[0, il] = (layer.Thickness / application_depth) * f_orgmat * amount
            elif(application_depth >= zmin and application_depth <= zmax):
                ORGMAT_am[0, il] = ((application_depth - zmin) / application_depth) * f_orgmat * amount
            elif(application_depth < zmin):
                ORGMAT_am[0, il] = 0
            else:
                pass
            CORG_am[0, il] = minip_C.calculate_organic_C(ORGMAT_am[0,il])

            if(cnratio == 0):
                NORG_am[0, il] = 0.
            else:
                NORG_am[0, il] = CORG_am[0, il] / cnratio
            zmin = zmax

        # Recreate the state variable matrices that contain both the existing amendment and the new amendments
        AGE = np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1]))
        AGE0 = np.zeros_like(AGE)
        ORGMAT = np.zeros_like(AGE)
        CORG = np.zeros_like(AGE)
        NORG = np.zeros_like(AGE)

        for am in range(0, AGE.shape[0]):
            for il in range(0, AGE.shape[1]):
                if(am < AGE.shape[0] - 1):
                    # Refill the matrix with the existing amendments
                    AGE0[am, il] = s.AGE0[am, il]
                    AGE[am, il] = s.AGE[am, il]
                    ORGMAT[am, il] = s.ORGMAT[am, il]
                    CORG[am, il] = s.CORG[am, il]
                    NORG[am, il] = s.NORG[am, il]
                else:
                    # Add the newly added amendment
                    AGE0[am, il] = AGE_am[0, il]
                    AGE[am, il] = AGE_am[0, il]
                    ORGMAT[am, il] = ORGMAT_am[0, il]
                    CORG[am, il] = CORG_am[0, il]
                    NORG[am, il] = NORG_am[0, il]

        s.AGE0 = AGE0
        s.AGE = AGE
        s.ORGMAT = ORGMAT
        s.CORG = CORG
        s.NORG = NORG

        # Add newly added ammonium and nitrate
        NH4 = np.zeros(len(self.soiln_profile))
        NO3 = np.zeros(len(self.soiln_profile))

        zmin = 0.
        for il in range(0, len(NH4)):
            layer = self.soiln_profile[il]
            if(application_depth > zmax):
                zmax = zmin + layer.Thickness
                NH4[il] = s.NH4[il] + (layer.Thickness / application_depth) * f_NH4N *  amount
                NO3[il] = s.NO3[il] + (layer.Thickness / application_depth) * f_NO3 *  amount
            elif(application_depth >= zmin and application_depth <= zmax):
                NH4[il] = s.NH4[il] + ((application_depth - zmin) / application_depth) * f_NH4N  * amount
                NO3[il] = s.NO3[il] + ((application_depth - zmin) / application_depth) * f_NO3  * amount
            elif(application_depth < zmin):
                NH4[il] = s.NH4[il]
                NO3[il] = s.NO3[il]
            else:
                pass
            zmin = zmax


        s.NH4= NH4
        s.NO3= NO3

        #r.unlock()
        ##r.FERT_N_SUPPLY = N_amount * N_recovery
        #r.lock()

    class SoilAmmoniumNModel():

        def calculate_mineralization_rate(self, dz, rNMINs_layer):
            RNH4MIN = (- rNMINs_layer).sum() / dz
            return RNH4MIN

        def calculate_nitrification_rate(self, cNH4, KNIT_REF, SM, SM0, T):
            fWNIT = self.calculate_soil_moisture_response_nitrification_rate_constant(SM, SM0)
            fT = self.calculate_temperature_response_nitrification_rate_constant(T)
            RNH4NIT = fWNIT * fT *  KNIT_REF * SM * cNH4
            return RNH4NIT

        def calculate_available_NH4(self, KSORP, NH4, RD, RHOD_kg_per_m2, SM, zmax, zmin):
            dz = zmax - zmin
            if(RD <= zmin):
                NH4_avail = 0.
            elif(RD > zmax):
                NH4_avail = (SM / ( KSORP * RHOD_kg_per_m2 + SM)) * NH4
            else:
                NH4_avail = ((RD - zmin)/ dz) * (SM / ( KSORP * RHOD_kg_per_m2 + SM)) * NH4
            return NH4_avail

        def calculate_NH4_plant_uptake_rate(self, KSORP, N_demand_soil, NH4, RD, RHOD_kg_per_m3, SM, zmax, zmin):
            return 0.

        def calculate_NH4_inflow_rate(self):
            return 0.

        def calculate_NH4_outflow_rate(self):
            return 0.

        def calculate_NH4_netflow(self):
            return 0.

        def calculate_soil_moisture_response_nitrification_rate_constant(self, SM, SM0):
            WFPS = SM / SM0
            fWNIT = 0.9 / (1. + np.exp(-15 *(WFPS - 0.45))) + 0.1 - 1/(1+np.exp(-50. * (WFPS - 0.95)))
            return fWNIT

        def calculate_temperature_response_nitrification_rate_constant(self, T):
            fT = 1/(1+np.exp(-0.26*(T-17.)))-1/(1+np.exp(-0.77*(T-41.9)))
            return fT

    class SoilNNitrateModel():
        def calculate_available_NO3(self, NO3, RD, zmax, zmin):
            dz = zmax - zmin
            if(RD <= zmin):
                NO3_avail_layer = 0.
            elif(RD > zmax):
                NO3_avail_layer = NO3            
            else:
                NO3_avail_layer = ((RD - zmin)/ dz) * NO3 
            return NO3_avail_layer

        def calculate_nitrification_rate(self, cNH4, KNIT_REF, SM, SM0):
            fWNIT = self.calculate_soil_moisture_response_nitrification_rate_constant(SM, SM0)
            RNH4NIT = fWNIT * KNIT_REF * SM * cNH4
            return RNH4NIT

        def calculate_denitrification_rate(self, cNO3, KDENIT_REF, MRCDIS, RCORGT, SM, SM0, T, WFPS_CRIT):
            fR = self.calculate_soil_respiration_response_denitrifiation_rate_constant(RCORGT, MRCDIS)
            fW = self.calculate_soil_moisture_response_denitrification_rate_constant(SM, SM0, WFPS_CRIT)
            fT = self.calculate_temperature_response_denitrification_rate_constant(T)
            RNO3DENIT = fW * fT * fR * KDENIT_REF * SM * cNO3
            return RNO3DENIT

        def calculate_NO3_plant_uptake_rate(self, N_demand_soil, NO3, RD, zmax, zmin):
            return 0.

        def calculate_NO3_inflow(self):
            return 0.

        def calculate_NO3_outflow(self):
            return 0.

        def calculate_NO3_netflow(self):
            return 0.

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
                rate = - k * OM
                return rate

            def calculate_dissimilation_rate_OM_T(self, OM, a, t, T):
                k = self.calculate_relative_dissimilation_rate_OM_T(a, t, T)
                rate = - k * OM
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
                rate = -k * Corg
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