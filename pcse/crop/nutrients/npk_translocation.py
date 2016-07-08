# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf

from ...traitlets import Float, Instance
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class NPK_Translocation(SimulationObject):
    """Does the bookkeeping for translocation of N/P/K from the roots, leaves
    and stems towards the storage organs of the crop.

    First the routine calculates the state of the translocatable amount of N/P/K.
    This translocatable amount is defined as the amount of N/P/K above the
    residual N/P/K amount calculated as the residual concentration times the
    living biomass. The residual amount is locked into the plant structural biomass
    and cannot be mobilized anymore. The translocatable amount is calculated for
    stems, roots and leaves and published as the state variables
    NTRANSLOCATABLE, PTRANSLOCATABLE and KTRANSLOCATABLE.

    The overal translocation rate is calculated as the minimum of supply (the
    translocatable amount) and demand from the storage organs as calculated in
    the component on Demand_Uptake.
    The actual rate of N/P/K translocation from the different plant organs is
    calculated assuming that the uptake rate is distributed over roots, stems and
    leaves in proportion to the translocatable amount for each organ.

    **Simulation parameters**

    ===============  ============================================= =======  ======================
     Name             Description                                   Type     Unit
    ===============  ============================================= =======  ======================
    NRESIDLV          Residual N fraction in leaves                 SCr     kg N kg-1 dry biomass
    PRESIDLV          Residual P fraction in leaves                 SCr     kg P kg-1 dry biomass
    KRESIDLV          Residual K fraction in leaves                 SCr     kg K kg-1 dry biomass

    NRESIDST          Residual N fraction in stems                  SCr     kg N kg-1 dry biomass
    PRESIDST          Residual P fraction in stems                  SCr     kg P kg-1 dry biomass
    KRESIDST          Residual K fraction in stems                  SCr     kg K kg-1 dry biomass

    NPK_TRANSLRT_FR   PK translocation from roots as a fraction     SCr     -
                      of resp. total NPK amounts translocated
                      from leaves and stems
    ===============  ============================================= ======= =======================


    **State variables**

    ================  ================================================= ===== ============
     Name               Description                                      Pbl      Unit
    ================  ================================================= ===== ============
    ATNLV              Translocatable N amount in living leaves           N    |kg N ha-1|
    ATPLV              Translocatable P amount in living leaves           N    |kg P ha-1|
    ATKLV              Translocatable K amount in living leaves           N    |kg K ha-1|
    ATNST              Translocatable N amount in living stems            N    |kg N ha-1|
    ATPST              Translocatable P amount in living stems            N    |kg P ha-1|
    ATKST              Translocatable K amount in living stems            N    |kg K ha-1|
    ATNRT              Translocatable N amount in living roots            N    |kg N ha-1|
    ATPRT              Translocatable P amount in living roots            N    |kg P ha-1|
    ATKRT              Translocatable K amount in living roots            N    |kg K ha-1|
    NTRANSLOCATABLE    Total N amount that can be translocated to the     Y    [kg N ha-1]
                       storage organs
    PTRANSLOCATABLE    Total P amount that can be translocated to the     Y    [kg P ha-1]
                       storage organs
    KTRANSLOCATABLE    Total K amount that can be translocated to the     Y    [kg K ha-1]
                       storage organs
    ================  ================================================= ===== ============


    **Rate variables**

    =======  ================================================= ==== ==============
     Name     Description                                      Pbl      Unit
    =======  ================================================= ==== ==============
    RNTLV     Weight increase (N) in leaves                     Y    |kg ha-1 d-1|
    RPTLV     Weight increase (P) in leaves                     Y    |kg ha-1 d-1|
    RKTLV     Weight increase (K) in leaves                     Y    |kg ha-1 d-1|
    RNTST     Weight increase (N) in stems                      Y    |kg ha-1 d-1|
    RPTST     Weight increase (P) in stems                      Y    |kg ha-1 d-1|
    RKTST     Weight increase (K) in stems                      Y    |kg ha-1 d-1|
    RNTRT     Weight increase (N) in roots                      Y    |kg ha-1 d-1|
    RPTRT     Weight increase (P) in roots                      Y    |kg ha-1 d-1|
    RKTRT     Weight increase (K) in roots                      Y    |kg ha-1 d-1|
    =======  ================================================= ==== ==============

    **Signals send or handled**

    None


    **External dependencies:**

    =========  =================================== =====================  ===========
     Name       Description                         Provided by            Unit
    =========  =================================== =====================  ===========
    DVS         Crop development stage              DVS_Phenology           -
    WST         Dry weight of living stems          WOFOST_Stem_Dynamics   |kg ha-1|
    WLV         Dry weight of living leaves         WOFOST_Leaf_Dynamics   |kg ha-1|
    WRT         Dry weight of living roots          WOFOST_Root_Dynamics   |kg ha-1|
    ANLV        Amount of N in leaves               NPK_Crop_Dynamics      |kg ha-1|
    ANST        Amount of N in stems                NPK_Crop_Dynamics      |kg ha-1|
    ANRT        Amount of N in roots                NPK_Crop_Dynamics      |kg ha-1|
    APLV        Amount of P in leaves               NPK_Crop_Dynamics      |kg ha-1|
    APST        Amount of P in stems                NPK_Crop_Dynamics      |kg ha-1|
    APRT        Amount of P in roots                NPK_Crop_Dynamics      |kg ha-1|
    AKLV        Amount of K in leaves               NPK_Crop_Dynamics      |kg ha-1|
    AKST        Amount of K in stems                NPK_Crop_Dynamics      |kg ha-1|
    AKRT        Amount of K in roots                NPK_Crop_Dynamics      |kg ha-1|
    =========  =================================== =====================  ===========
    """

    class Parameters(ParamTemplate):
        NRESIDLV = Float(-99.)  # residual N fraction in leaves [kg N kg-1 dry biomass]
        NRESIDST = Float(-99.)  # residual N fraction in stems [kg N kg-1 dry biomass]
        NRESIDRT = Float(-99.)  # residual N fraction in roots [kg N kg-1 dry biomass]

        PRESIDLV = Float(-99.)  # residual P fraction in leaves [kg P kg-1 dry biomass]
        PRESIDST = Float(-99.)  # residual P fraction in stems [kg P kg-1 dry biomass]
        PRESIDRT = Float(-99.)  # residual P fraction in roots [kg P kg-1 dry biomass]

        KRESIDLV = Float(-99.)  # residual K fraction in leaves [kg P kg-1 dry biomass]
        KRESIDST = Float(-99.)  # residual K fraction in stems [kg P kg-1 dry biomass]
        KRESIDRT = Float(-99.)  # residual K fraction in roots [kg P kg-1 dry biomass]

        NPK_TRANSLRT_FR = Float(-99.)  # NPK translocation from roots as a fraction of
                                       # resp. total NPK amounts translocated from leaves
                                       # and stems

    class RateVariables(RatesTemplate):
        RNTLV = Float(-99.)  # N translocation rate from leaves [kg ha-1 d-1]
        RNTST = Float(-99.)  # N translocation rate from stems [kg ha-1 d-1]
        RNTRT = Float(-99.)  # N translocation rate from roots [kg ha-1 d-1]

        RPTLV = Float(-99.)  # P translocation rate from leaves [kg ha-1 d-1]
        RPTST = Float(-99.)  # P translocation rate from stems [kg ha-1 d-1]
        RPTRT = Float(-99.)  # P translocation rate from roots [kg ha-1 d-1]

        RKTLV = Float(-99.)  # K translocation rate from leaves [kg ha-1 d-1]
        RKTST = Float(-99.)  # K translocation rate from stems [kg ha-1 d-1]
        RKTRT = Float(-99.)  # K translocation rate from roots [kg ha-1 d-1]

    class StateVariables(StatesTemplate):
        ATNLV = Float(-99.)  # translocatable N amount in leaves [kg N ha-1]
        ATNST = Float(-99.)  # translocatable N amount in stems [kg N ha-1]
        ATNRT = Float(-99.)  # translocatable N amount in roots [kg N ha-1]
        
        ATPLV = Float(-99.)  # translocatable P amount in leaves [kg N ha-1]
        ATPST = Float(-99.)  # translocatable P amount in stems [kg N ha-1]
        ATPRT = Float(-99.)  # translocatable P amount in roots [kg N ha-1]
        
        ATKLV = Float(-99.)  # translocatable K amount in leaves [kg N ha-1
        ATKST = Float(-99.)  # translocatable K amount in stems [kg N ha-1]
        ATKRT = Float(-99.)  # translocatable K amount in roots [kg N ha-1]

        NTRANSLOCATABLE = Float(-99.)  # Total N amount that can be translocated to the storage organs [kg N ha-1]
        PTRANSLOCATABLE = Float(-99.)  # Total P amount that can be translocated to the storage organs [kg P ha-1]
        KTRANSLOCATABLE = Float(-99.)  # Total K amount that can be translocated to the storage organs [kg K ha-1]

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE instance
        :param parvalues: dictionary with WOFOST cropdata key/value pairs
        """

        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(kiosk, publish=["RNTLV", "RNTST", "RNTRT",
                                                        "RPTLV", "RPTST", "RPTRT",
                                                        "RKTLV", "RKTST", "RKTRT"])

        self.states = self.StateVariables(kiosk,
            ATNLV=0., ATNST=0., ATNRT=0., ATPLV=0., ATPST=0., ATPRT=0., ATKLV=0., ATKST=0. ,ATKRT=0.,
            NTRANSLOCATABLE=0., PTRANSLOCATABLE=0., KTRANSLOCATABLE=0.,
            publish=["NTRANSLOCATABLE", "PTRANSLOCATABLE", "KTRANSLOCATABLE"])
        self.kiosk = kiosk
        
    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states

        RNUSO = self.kiosk["RNUSO"]  # N uptake storage organs
        RPUSO = self.kiosk["RPUSO"]  # P uptake storage organs
        RKUSO = self.kiosk["RKUSO"]  # K uptake storage organs

#       partionioning of the uptake
        # if max amount of translocatable N = 0 then
        # N translocation rate is 0
        if s.NTRANSLOCATABLE > 0.:
            r.RNTLV = RNUSO * s.ATNLV / s.NTRANSLOCATABLE
            r.RNTST = RNUSO * s.ATNST / s.NTRANSLOCATABLE
            r.RNTRT = RNUSO * s.ATNRT / s.NTRANSLOCATABLE
        else:
            r.RNTLV = r.RNTST = r.RNTRT = 0.

        # if max amount of translocatable P = 0 then
        # P translocation rate is 0
        if s.PTRANSLOCATABLE > 0:
            r.RPTLV = RPUSO * s.ATPLV / s.PTRANSLOCATABLE
            r.RPTST = RPUSO * s.ATPST / s.PTRANSLOCATABLE
            r.RPTRT = RPUSO * s.ATPRT / s.PTRANSLOCATABLE
        else:
            r.RPTLV = r.RPTST = r.RPTRT = 0.

        # if max amount of translocatable K = 0 then
        # K translocation rate is 0
        if s.KTRANSLOCATABLE > 0:
            r.RKTLV = RKUSO * s.ATKLV / s.KTRANSLOCATABLE
            r.RKTST = RKUSO * s.ATKST / s.KTRANSLOCATABLE
            r.RKTRT = RKUSO * s.ATKRT / s.KTRANSLOCATABLE
        else:
            r.RKTLV = r.RKTST = r.RKTRT = 0.

    @prepare_states
    def integrate(self, day, delt=1.0):
        p = self.params
        s = self.states
        
        WLV = self.kiosk["WLV"]
        WST = self.kiosk["WST"]
        WRT = self.kiosk["WRT"]
        
        ANLV = self.kiosk["ANLV"]
        ANST = self.kiosk["ANST"]
        ANRT = self.kiosk["ANRT"]
        
        APLV = self.kiosk["APLV"]
        APST = self.kiosk["APST"]
        APRT = self.kiosk["APRT"]
        
        AKLV = self.kiosk["AKLV"]
        AKST = self.kiosk["AKST"]
        AKRT = self.kiosk["AKRT"]

#       translocatable N amount in the organs [kg N ha-1]
        s.ATNLV = max(0., ANLV - WLV * p.NRESIDLV)
        s.ATNST = max(0., ANST - WST * p.NRESIDST)
        s.ATNRT = max((s.ATNLV + s.ATNST) * p.NPK_TRANSLRT_FR, ANRT - WRT * p.NRESIDRT)

#       translocatable P amount in the organs [kg P ha-1]
        s.ATPLV = max(0., APLV - WLV * p.PRESIDLV)
        s.ATPST = max(0., APST - WST * p.PRESIDST)
        s.ATPRT = max((s.ATPLV + s.ATPST) * p.NPK_TRANSLRT_FR, APRT - WRT * p.PRESIDRT)

#       translocatable K amount in the organs [kg K ha-1]
        s.ATKLV = max(0., AKLV - WLV * p.KRESIDLV)
        s.ATKST = max(0., AKST - WST * p.KRESIDST)
        s.ATKRT = max((s.ATKLV + s.ATKST) * p.NPK_TRANSLRT_FR, AKRT - WRT * p.KRESIDRT)

#       total translocatable NPK amount in the organs [kg N ha-1]
        s.NTRANSLOCATABLE = s.ATNLV + s.ATNST + s.ATNRT
        s.PTRANSLOCATABLE = s.ATPLV + s.ATPST + s.ATPRT
        s.KTRANSLOCATABLE = s.ATKLV + s.ATKST + s.ATKRT
