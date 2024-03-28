# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
from collections import namedtuple

from ...base import StatesTemplate, ParamTemplate, SimulationObject, RatesTemplate
from ...decorators import prepare_rates, prepare_states
from ...traitlets import HasTraits, Float, Int, Instance
from ...util import AfgenTrait

MaxNutrientConcentrations = namedtuple("MaxNutrientConcentrations",
                                       ["NMAXLV","NMAXST", "NMAXRT", "NMAXSO"])

class N_Demand_Uptake(SimulationObject):
    """Calculates the crop N demand and its uptake from the soil.

    Crop N demand is calculated as the difference between the actual N
    (kg N per kg biomass) in the vegetative plant organs (leaves, stems and roots)
    and the maximum N concentration for each organ. N uptake is then 
    estimated as the minimum of supply from the soil and demand from the crop.

    Nitrogen fixation (leguminous plants) is calculated by assuming that a
    fixed fraction of the daily N demand is supplied by nitrogen fixation.
    The remaining part has to be supplied by the soil.

    The N demand of the storage organs is calculated in a somewhat
    different way because it is assumed that the demand from the storage
    organs is fulfilled by translocation of N/P/K from the leaves, stems
    and roots. Therefore the uptake of the storage organs is calculated
    as the minimum of the daily translocatable N supply and the demand from
    the storage organs.

    **Simulation parameters**

    ============  =============================================  ======================
     Name          Description                                    Unit
    ============  =============================================  ======================
    NMAXLV_TB      Maximum N concentration in leaves as          kg N kg-1 dry biomass
                   function of DVS
    NMAXRT_FR      Maximum N concentration in roots as fraction  -
                   of maximum N concentration in leaves
    NMAXST_FR      Maximum N concentration in stems as fraction  -
                   of maximum N concentration in leaves
    NMAXSO         Maximum N concentration in storage organs     kg N kg-1 dry biomass
    NCRIT_FR       Critical N concentration as fraction of       -
                   maximum N concentration for vegetative
                   plant organs as a whole (leaves + stems)
    TCNT           Time coefficient for N translation to         days
                   storage organs
    NFIX_FR        fraction of crop nitrogen uptake by           kg N kg-1 dry biomass biological fixation
    RNUPTAKEMAX    Maximum rate of N uptake                      |kg N ha-1 d-1|
    ============  =============================================  ======================

    **State variables**

    ============= ================================================= ==== ============
     Name          Description                                      Pbl      Unit
    ============= ================================================= ==== ============
    NuptakeTotal  Total N uptake by the crop                        N   |kg N ha-1|
    NfixTotal     Total N fixated by the crop                       N   |kg N ha-1|
    NdemandST     N Demand in living stems                          N   |kg N ha-1|
    NdemandRT     N Demand in living roots                          N   |kg N ha-1|
    NdemandSO     N Demand in storage organs                        N   |kg N ha-1|
    ==========    ================================================= ==== ============


    **Rate variables**

    ===========  ================================================= ==== ================
     Name         Description                                      Pbl      Unit
    ===========  ================================================= ==== ================
    RNuptakeLV     Rate of N uptake in leaves                        Y   |kg N ha-1 d-1|
    RNuptakeST     Rate of N uptake in stems                         Y   |kg N ha-1 d-1|
    RNuptakeRT     Rate of N uptake in roots                         Y   |kg N ha-1 d-1|
    RNuptakeSO     Rate of N uptake in storage organs                Y   |kg N ha-1 d-1|
    RNuptake       Total rate of N uptake                            Y   |kg N ha-1 d-1|
    RNfixation     Rate of N fixation                                Y   |kg N ha-1 d-1|
    NdemandLV      N Demand in living leaves                         N   |kg N ha-1|
    NdemandST      N Demand in living stems                          N   |kg N ha-1|
    NdemandRT      N Demand in living roots                          N   |kg N ha-1|
    NdemandSO      N Demand in storage organs                        N   |kg N ha-1|
    Ndemand        Total crop N demand                               N   |kg N ha-1 d-1|
    ===========  ================================================= ==== ================

    **Signals send or handled**

    None

    **External dependencies**

    ================  =================================== ====================  ===========
     Name              Description                        Provided by            Unit
    ================  =================================== ====================  ===========
    DVS               Crop development stage              DVS_Phenology              -
    TRA               Crop transpiration                  Evapotranspiration     |cm d-1|
    TRAMX             Potential crop transpiration        Evapotranspiration     |cm d-1|
    NAVAIL            Total available N from soil         N_Soil_Dynamics      |kg ha-1|
    ================  =================================== ====================  ===========

    """

    class Parameters(ParamTemplate):
        NMAXLV_TB = AfgenTrait()  # maximum N concentration in leaves as function of dvs        
        DVS_N_TRANSL = Float(-99.)

        NMAXRT_FR = Float(-99.)  # maximum N concentration in roots as fraction of maximum N concentration in leaves
        NMAXST_FR = Float(-99.)  # maximum N concentration in stems as fraction of maximum N concentration in leaves        
        NMAXSO = Float(-99.)  # maximum P concentration in storage organs [kg N kg-1 dry biomass]        
        TCNT = Float(-99.)  # time coefficient for N translocation to storage organs [days]

        NFIX_FR = Float(-99.)  # fraction of crop nitrogen uptake by biological fixation
        RNUPTAKEMAX = Float()  # Maximum N uptake rate
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NRESIDRT = Float(-99.)  # residual N fraction in roots [kg N kg-1 dry biomass]

    class RateVariables(RatesTemplate):
        RNtranslocationLV = Float(-99.)  # N translocation rate from leaves [kg ha-1 d-1]
        RNtranslocationST = Float(-99.)  # N translocation rate from stems [kg ha-1 d-1]
        RNtranslocationRT = Float(-99.)  # N translocation rate from roots [kg ha-1 d-1]
        RNtranslocation = Float(-99.)    # N translocation rate to storage organs [kg ha-1 d-1]

        RNuptakeLV = Float(-99.)  # N uptake rates in organs [kg ha-1 d -1]
        RNuptakeST = Float(-99.)
        RNuptakeRT = Float(-99.)
        RNuptakeSO = Float(-99.)

        RNuptake = Float(-99.)  # Total N uptake rates [kg ha-1 d -1]
        RNfixation = Float(-99.)  # Total N fixated

        NdemandLV = Float(-99.)  # N demand in organs [kg ha-1]
        NdemandST = Float(-99.)
        NdemandRT = Float(-99.)
        NdemandSO = Float(-99.)

        Ndemand = Float()  # Total N/P/K demand of the crop

    class StateVariables(StatesTemplate):
        NtranslocatableLV = Float(-99.)  # translocatable N amount in leaves [kg N ha-1]
        NtranslocatableST = Float(-99.)  # translocatable N amount in stems [kg N ha-1]
        NtranslocatableRT = Float(-99.)  # translocatable N amount in roots [kg N ha-1]
        Ntranslocatable = Float(-99.)

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: a ParameterProvider with parameter key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.kiosk = kiosk

        self.rates = self.RateVariables(kiosk,
            publish=["RNtranslocationLV", "RNtranslocationST", "RNtranslocationRT", "RNtranslocation", 
                     "RNuptakeLV", "RNuptakeST", "RNuptakeRT", "RNuptakeSO","RNuptake", "RNfixation"])

        self.states = self.StateVariables(kiosk, NtranslocatableLV=0., NtranslocatableST=0., NtranslocatableRT=0., 
                                          Ntranslocatable=0., publish=["Ntranslocatable"])

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        p = self.params
        k = self.kiosk
        s = self.states

        delt = 1.0
        mc = self._compute_N_max_concentrations()

        # No nutrients are absorbed when severe water shortage occurs i.e. RFTRA <= 0.01
        if k.RFTRA > 0.01:
            NutrientLIMIT = 1.0
        else:
            NutrientLIMIT = 0.

        # N demand [kg ha-1]
        r.NdemandLV = max(mc.NMAXLV * k.WLV - k.NamountLV, 0.) + max(k.GRLV * mc.NMAXLV, 0) * delt
        r.NdemandST = max(mc.NMAXST * k.WST - k.NamountST, 0.) + max(k.GRST * mc.NMAXST, 0) * delt
        r.NdemandRT = max(mc.NMAXRT * k.WRT - k.NamountRT, 0.) + max(k.GRRT * mc.NMAXRT, 0) * delt
        r.NdemandSO = max(mc.NMAXSO * k.WSO - k.NamountSO, 0.) + max(k.GRSO * mc.NMAXSO, 0) * delt

        r.Ndemand = r.NdemandLV + r.NdemandST + r.NdemandRT + r.NdemandSO


        # biological nitrogen fixation
        r.RNfixation = (max(0., p.NFIX_FR * r.Ndemand) * NutrientLIMIT)

        # Calculate translocatable nitrogen in different organs
        if(k.DVS < p.DVS_N_TRANSL):
            s.NTranslocatableLV = 0.
            s.NTranslocatableRT = 0.
            s.NTranslocatableST = 0.
            s.NTranslocatable = 0.
        else:
            s.NTranslocatableLV = max(0., k.NamountLV - k.WLV * p.NRESIDLV)
            s.NTranslocatableRT = max(0., k.NamountRT - k.WRT * p.NRESIDRT)
            s.NTranslocatableST = max(0., k.NamountST - k.WST * p.NRESIDST)
            s.NTranslocatable = s.NTranslocatableLV + s.NTranslocatableRT + s.NTranslocatableST

        r.RNtranslocation = min(r.NdemandSO/delt, s.NTranslocatable / p.TCNT)

        if(s.NTranslocatable == 0):
            r.RNtranslocationLV = 0.
            r.RNtranslocationRT = 0.
            r.RNtranslocationST = 0.
        else:
            r.RNtranslocationLV = r.RNtranslocation * (s.NTranslocatableLV / s.NTranslocatable)
            r.RNtranslocationRT = r.RNtranslocation * (s.NTranslocatableRT / s.NTranslocatable) 
            r.RNtranslocationST = r.RNtranslocation * (s.NTranslocatableST / s.NTranslocatable)

        r.RNuptake = (max(0., min(r.Ndemand - r.RNfixation, k.NAVAIL, p.RNUPTAKEMAX)) * NutrientLIMIT)

        if r.Ndemand == 0:
            r.RNuptakeLV = 0.
            r.RNuptakeRT = 0.
            r.RNuptakeST = 0.
            r.RNuptakeSO = 0.
        else:
            r.RNuptakeLV = max(0.,min(r.NdemandLV/delt + r.RNtranslocationLV, r.RNuptake * (r.NdemandLV/delt + r.RNtranslocationLV) / r.Ndemand))
            r.RNuptakeRT = max(0.,min(r.NdemandRT/delt + r.RNtranslocationRT, r.RNuptake * (r.NdemandRT/delt + r.RNtranslocationRT) / r.Ndemand))
            r.RNuptakeST = max(0.,min(r.NdemandST/delt + r.RNtranslocationST, r.RNuptake * (r.NdemandST/delt + r.RNtranslocationST) / r.Ndemand))
            r.RNuptakeSO = max(0.,min(r.NdemandSO/delt - r.RNtranslocation,   r.RNuptake * (r.NdemandSO/delt - r.RNtranslocation) / r.Ndemand))


    @prepare_states
    def integrate(self, day, delt=1.0):
        pass

    def _compute_N_max_concentrations(self):
        """Computes the maximum N concentrations in leaves, stems, roots and storage organs.
        
        Note that max concentrations are first derived from the dilution curve for leaves. 
        Maximum concentrations for stems and roots are computed as a fraction of the 
        concentration for leaves.
        """

        p = self.params
        k = self.kiosk
        NMAXLV = p.NMAXLV_TB(k.DVS)

        max_N_conc = MaxNutrientConcentrations(
            # Maximum N concentrations in leaves [kg N kg-1 DM]
            NMAXLV=NMAXLV,
            # Maximum N concentrations in stems and roots [kg N kg-1 DM]
            NMAXST=(p.NMAXST_FR * NMAXLV),
            NMAXRT=p.NMAXRT_FR * NMAXLV,
            NMAXSO=p.NMAXSO
        )

        return max_N_conc
