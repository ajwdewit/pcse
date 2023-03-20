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
        AGE0   = Instance(np.matrix) # Initial age of material (d)
        AGE    = Instance(np.matrix) # Appearant age of material (d)
        ORGMAT = Instance(np.matrix) # Amount of organic matter (kg ORG ha-1)
        CORG   = Instance(np.matrix) # Amount of C in organic matter (kg C ha-1)
        NORG   = Instance(np.matrix) # 
        NH4    = Instance(np.ndarray) # Amount of NH4-N (kg N ha-1)
        NO3    = Instance(np.ndarray) # Amount of NO3-N (kg N ha-1)
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1

        ORGMATT = Float()
        CORGT = Float()       
        NORGT = Float()
        RMINT = Float()
        NH4T = Float()
        NO3T = Float

    class RateVariables(RatesTemplate):
        RAGE = Instance(np.matrix)
        RORGMAT = Instance(np.matrix)
        RCORG = Instance(np.matrix)
        RNORG = Instance(np.matrix)
        
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
        AGE =  np.matrix(np.zeros(len(self.soiln_profile)))
        AGE0 = np.matrix(np.zeros(len(self.soiln_profile)))
        ORGMAT = np.matrix(np.zeros(len(self.soiln_profile)))
        CORG =  np.matrix(np.zeros(len(self.soiln_profile)))
        NORG =  np.matrix(np.zeros(len(self.soiln_profile)))
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            NH4[il] = 0.
            NO3[il] = 0.
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

        states = {"NAVAIL": 100., "NO3": NO3, "NH4": NH4, "AGE": AGE, "AGE0": AGE0, "ORGMAT": ORGMAT, "CORG": CORG, "NORG": NORG,  
                  "ORGMATT": ORGMATT,  "CORGT": CORGT, "NORGT": NORGT, "RMINT": RMINT, "NH4T": NH4T, "NO3T": NO3T}
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
        r.RAGE = np.matrix(np.zeros((s.AGE.shape[0],s.AGE.shape[1])))
        r.RAGE = np.matrix(np.zeros((s.AGE.shape[0],s.AGE.shape[1])))
        r.RORGMAT = np.matrix(np.zeros((s.AGE.shape[0],s.AGE.shape[1])))
        r.RCORG = np.matrix(np.zeros((s.AGE.shape[0],s.AGE.shape[1])))
        r.RNORG = np.matrix(np.zeros((s.AGE.shape[0],s.AGE.shape[1])))

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
        r.RNH4MIN = np.zeros(len(self.soiln_profile))
        r.RNH4NITR = np.zeros(len(self.soiln_profile))
        r.RNH4UP = np.zeros(len(self.soiln_profile))
        r.RNH4IN = np.zeros(len(self.soiln_profile))
        r.RNH4OUT = np.zeros(len(self.soiln_profile))

        samm = self.SoilAmmoniumNModel()
        for il in range(0, len(self.soiln_profile)):
            dz = self.soiln_profile[il].Thickness * self.cm_to_m
            SM0 = self.soiln_profile[il].SM0
            RNMIN_kg_per_m3 = r.RNORG[:,il] * self.m2_to_ha
            cNH4 = s.NH4[il] * self.m2_to_ha / (dz * k.SM[il])

            r.RNH4MIN[il] = samm.calculate_mineralization_rate(dz, RNMIN_kg_per_m3)
            r.RNH4NITR[il] = samm.calculate_nitrification_rate(cNH4, p.KNIT_REF, k.SM[il],SM0, T)
            r.RNH4[il] = (1/self.m2_to_ha) * dz * (r.RNH4MIN[il] - r.RNH4NITR[il] - r.RNH4UP[il])

        # Initialize rates for nitrate
        r.RNO3  = np.zeros(len(self.soiln_profile))
        r.RNO3NITR = np.zeros(len(self.soiln_profile))
        r.RNO3DENITR = np.zeros(len(self.soiln_profile))
        r.RNO3UP = np.zeros(len(self.soiln_profile))
        r.RNO3IN = np.zeros(len(self.soiln_profile))
        r.RNO3OUT = np.zeros(len(self.soiln_profile))

        sni = self.SoilNNitrateModel()
        for il in range(0, len(self.soiln_profile)):
            dz = self.soiln_profile[il].Thickness * self.cm_to_m
            cNO3 = s.NO3[il] * self.m2_to_ha / (dz * k.SM[il])
            RCORGT_kg_per_m2 = - r.RCORG.sum() * self.ha_to_m2
            SM0 = self.soiln_profile[il].SM0
            r.RNO3NITR[il] = r.RNH4NITR[il]
            r.RNO3DENITR[il] = sni.calculate_denitrification_rate(cNO3, p.KDENIT_REF, p.MRCDIS, RCORGT_kg_per_m2, k.SM[il], SM0, T, p.WFPS_CRIT)
            r.RNO3[il] =  (1/self.m2_to_ha) * dz *  (r.RNO3NITR[il]  - r.RNO3DENITR[il])

    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        s = self.states

        s.NAVAIL = s.NAVAIL

        AGE =  np.matrix(np.zeros((s.AGE.shape[0], s.AGE.shape[1])))
        ORGMAT = np.matrix(np.zeros((s.AGE.shape[0], s.AGE.shape[1])))
        CORG = np.matrix(np.zeros((s.AGE.shape[0], s.AGE.shape[1])))
        NORG = np.matrix(np.zeros((s.AGE.shape[0], s.AGE.shape[1])))
        NO3 = np.zeros(len(self.soiln_profile))
        NH4 = np.zeros(len(self.soiln_profile))

        for am in range(0, r.RAGE.shape[0]):
            for il in range(0, r.RAGE.shape[1]):
                AGE[am, il] = s.AGE[am,il] + r.RAGE[am, il] * delt
                ORGMAT[am, il] = s.ORGMAT[am,il] + r.RORGMAT[am, il] * delt
                CORG[am, il] = s.CORG[am,il] + r.RCORG[am, il] * delt
                NORG[am, il] = s.NORG[am, il] + r.RNORG[am, il] * delt

        for il in range(0, len(s.NH4)):
            NH4[il] = s.NH4[il] + r.RNH4[il] * delt
            NO3[il] = s.NO3[il] + r.RNO3[il] * delt

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


    def _on_APPLY_N(self, amount=None, application_depth = None, cnratio=None, f_orgmat=None, f_NH4N = None, f_NO3 = None, initial_age =None):
        r = self.rates
        s = self.states

        AGE0_am =  np.matrix(np.zeros(len(self.soiln_profile)))
        AGE_am = np.matrix(np.zeros(len(self.soiln_profile)))
        ORGMAT_am = np.matrix(np.zeros(len(self.soiln_profile)))
        CORG_am =  np.matrix(np.zeros(len(self.soiln_profile)))
        NORG_am =  np.matrix(np.zeros(len(self.soiln_profile)))

        # Make matrix for new amendment
        zmin = 0.
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            AGE0_am[0,il] = initial_age * self.y_to_d
            AGE_am[0,il] = initial_age * self.y_to_d            
            zmax = zmin + self.soiln_profile[il].Thickness
            if(application_depth > zmax):
                ORGMAT_am[0, il] = (layer.Thickness / application_depth) * amount
            elif(application_depth >= zmin and application_depth <= zmax):
                ORGMAT_am[0, il] = ((application_depth - zmin) / application_depth) * amount
            elif(application_depth < zmin):
                ORGMAT_am[0, il] = 0
            else:
                pass
            CORG_am[0, il] = minip_C.calculate_organic_C(ORGMAT_am[0,il])
            NORG_am[0, il] = CORG_am[0, il] / cnratio
            zmin = zmax

        # Recreate the state variable matrices that contain both the existing amendment and the new amendments
        AGE0 = np.matrix(np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1])))
        AGE = np.matrix(np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1])))
        ORGMAT = np.matrix(np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1])))
        CORG = np.matrix(np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1])))
        NORG = np.matrix(np.zeros((s.AGE.shape[0] + 1,s.AGE.shape[1])))

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

        ## Add newly added ammonium and nitrate
        #NH4_amm = np.array(len(self.soiln_profile))
        #NO3_amm = np.array(len(self.soiln_profile))

        #for il in range(0, len(NH4_amm)):
        #    if(application_depth > zmax):
        #        NH4_amm[il] = (layer.Thickness / application_depth) * amount
        #        NO3_amm[il] = (layer.Thickness / application_depth) * amount
        #    elif(application_depth >= zmin and application_depth <= zmax):
        #        NH4_amm[il] = ((application_depth - zmin) / application_depth) * amount
        #        NO3_amm[il] = (layer.Thickness / application_depth) * amount
        #    elif(application_depth < zmin):
        #        NH4_amm[il] = 0
        #        NO3_amm[il] = (layer.Thickness / application_depth) * amount
        #    else:
        #        pass

        #s.NH4+= NH4_amm
        #s.NO3+= NO3_amm

        r.unlock()
        #r.FERT_N_SUPPLY = N_amount * N_recovery
        r.lock()

    class SoilAmmoniumNModel():

        def calculate_mineralization_rate(self, dz, rNMINs_layer):
            RNH4MIN = (- rNMINs_layer).sum() / dz
            return RNH4MIN

        def calculate_nitrification_rate(self, cNH4, KNIT_REF, SM, SM0, T):
            fWNIT = self.calculate_soil_moisture_response_nitrification_rate_constant(SM, SM0)
            fT = self.calculate_temperature_response_nitrification_rate_constant(T)
            RNH4NIT = fWNIT * fT *  KNIT_REF * SM * cNH4
            return RNH4NIT

        def calculate_NH4_plant_uptake_rate(self):
            return 0.

        def calculate_NH4_inflow_rate(self):
            return 0.

        def calculate_NH4_outflow_rate(self):
            return 0.

        def calculate_NH4_netflow(self):
            return 0.

        def calculate_soil_moisture_response_nitrification_rate_constant(self, SM, SM0):
            WFPS = SM / SM0
            fWNIT = 0.9 / (1 + np.exp(-15 *(WFPS - 0.45))) + 0.1 - 1/(1+np.exp(-50 * (WFPS - 0.95)))
            return fWNIT

        def calculate_temperature_response_nitrification_rate_constant(self, T):
            fT = 1/(1+np.exp(-0.26*(T-17)))-1/(1+np.exp(-0.77*(T-41.9)))
            return fT

    class SoilNNitrateModel():
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

        def calculate_NO3_plant_uptake_rate(self):
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