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
    cm2_to_ha = 1e-8
    y_to_d = 365.25

    # placeholder for soil object
    soiln_profile = None

    class StateVariables(StatesTemplate):
        AGE0   = Instance(np.matrix) # Initial age of material (d)
        AGE    = Instance(np.matrix) # Appearant age of material (d)
        ORGMAT = Instance(np.matrix) # Amount of organic matter (kg ORG ha-1)
        CORG   = Instance(np.matrix) # Amount of C in organic matter (kg C ha-1)
        NORG   = Instance(np.matrix) # 
        NAVAIL = Float(-99.)  # total mineral N from soil and fertiliser  kg N ha-1

        NORGT = Float()

    class RateVariables(RatesTemplate):
        RAGE = Instance(np.matrix)
        RORGMAT = Instance(np.matrix)
        RCORG = Instance(np.matrix)
        RNORG = Instance(np.matrix)        

    class Parameters(ParamTemplate):
        A0SOM = Float()            # Initial age of humus
        CNRatioBio = Float()       # C:N Ratio of microbial biomass
        FASDIS = Float()           # Fraction of assimilation to dissimilation

    def initialize(self, day, kiosk, parvalues):
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)
        self.soiln_profile = SoilNProfile(parvalues)

        NORGT = 0.

        # Initialize states
        parvalues._soildata["soil_profile"] = self.soiln_profile
        AGE =  np.matrix(np.zeros(len(self.soiln_profile)))
        AGE0 = np.matrix(np.zeros(len(self.soiln_profile)))
        ORGMAT = np.matrix(np.zeros(len(self.soiln_profile)))
        CORG =  np.matrix(np.zeros(len(self.soiln_profile)))
        NORG =  np.matrix(np.zeros(len(self.soiln_profile)))
        minip_C = self.SoilOrganicNModel.MINIP_C()
        for il, layer in enumerate(self.soiln_profile):
            AGE0[0,il] = self.params.A0SOM * self.y_to_d
            AGE[0,il] = self.params.A0SOM * self.y_to_d
            ORGMAT[0,il] = layer.RHOD * layer.FSOMI * layer.Thickness *  self.g_to_kg / self.cm2_to_ha
            CORG[0,il] = minip_C.calculate_organic_C(ORGMAT[0,il])
            NORG[0,il] = CORG[0, il] / layer.CNRatioSOMI
        states = {"NAVAIL": 100., "AGE": AGE, "AGE0": AGE0, "ORGMAT": ORGMAT, "CORG": CORG, "NORG": NORG, "NORGT": NORGT}
        self.states = self.StateVariables(kiosk, publish=["NAVAIL", "ORGMAT", "CORG", "NORG"], **states)
        self.rates = self.RateVariables(kiosk)

        self._connect_signal(self._on_APPLY_N, signals.apply_n)

    @prepare_rates
    def calc_rates(self, day, drv):
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
                r.RORGMAT[am,il] = janssen.calculate_dissimilation_rate_OM(s.ORGMAT[am,il], s.AGE0[am,il], s.AGE[am,il])
                r.RCORG[am,il] = minip_c.calculate_dissimilation_rate_C(janssen, s.ORGMAT[am,il], s.AGE0[am,il], s.AGE[am,il], T)
                r.RNORG[am,il] = minip_n.calculate_dissimilation_rate_N(janssen, minip_c, s.ORGMAT[am,il], s.NORG[am,il], s.AGE0[am,il], p.FASDIS, p.CNRatioBio, s.AGE[am,il], T)

    @prepare_states
    def integrate(self, day, delt=1.0):
        r = self.rates
        s = self.states

        for am in range(0, r.RAGE.shape[0]):
            for il in range(0, r.RAGE.shape[1]):                
                self.states.ORGMAT[am, il] += r.RORGMAT[am, il] * delt
                self.states.CORG[am, il] += r.RCORG[am, il] * delt
                self.states.NORG[am, il] += r.RNORG[am, il] * delt
                self.touch()

        self.states.NORGT += -np.sum(r.RNORG) * delt

    def _on_APPLY_N(self, amount=None, application_depth = None, cnratio=None, f_orgmat=None, f_NH4N = None, f_NO3 = None, initial_age =None):
        r = self.rates
        s = self.states
        r.unlock()
        #r.FERT_N_SUPPLY = N_amount * N_recovery
        r.lock()

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