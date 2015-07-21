# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Allard de Wit and Iwan Supit (allard.dewit@wur.nl), July 2015
# Approach based on LINTUL N/P/K made by Joost Wolf

from ...traitlets import Float, Instance
from ...decorators import prepare_rates, prepare_states
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject

class NPK_Translocation(SimulationObject):
    """Calculates the translocatable amounts of N/P/K for the different organs.

    The translocatable amount of N/P/K is defined as the amount above the
    residual concentration which is locked into the plant structural biomass
    and which cannot be mobilized anymore.

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

    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk translocation __call__()
        """

        self.params = self.Parameters(cropdata)
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

#       max amount translocatable NPK [kg ha-1 d-1]
        ATN = s.ATNLV + s.ATNST + s.ATNRT
        ATP = s.ATPLV + s.ATPST + s.ATPRT
        ATK = s.ATKLV + s.ATKST + s.ATKRT

#       partionioning of the uptake
        # if max amount of translocatable N = 0 then
        # N translocation rate is 0
        if ATN > 0.:
            r.RNTLV = RNUSO * s.ATNLV / ATN
            r.RNTST = RNUSO * s.ATNST / ATN
            r.RNTRT = RNUSO * s.ATNRT / ATN
        else:
            r.RNTLV = r.RNTST = r.RNTRT = 0.

        # if max amount of translocatable P = 0 then
        # P translocation rate is 0
        if ATP > 0:
            r.RPTLV = RPUSO * s.ATPLV / ATP
            r.RPTST = RPUSO * s.ATPST / ATP
            r.RPTRT = RPUSO * s.ATPRT / ATP
        else:
            r.RPTLV = r.RPTST = r.RPTRT = 0.

        # if max amount of translocatable K = 0 then
        # K translocation rate is 0
        if ATK > 0:
            r.RKTLV = RKUSO * s.ATKLV / ATK
            r.RKTST = RKUSO * s.ATKST / ATK
            r.RKTRT = RKUSO * s.ATKRT / ATK
        else:
            r.RKTLV = r.RKTST = r.RKTRT = 0.

    @prepare_states
    def integrate(self, day):
        p = self.params
        s = self.states
        
        DVS = self.kiosk["DVS"]
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
