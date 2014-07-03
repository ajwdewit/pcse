#!/usr/bin/env python

from ...decorators import prepare_rates, prepare_states
from ...traitlets import Float, Instance, AfgenTrait
from ...base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject, VariableKiosk

class npk_losses(SimulationObject):
    
    class RateVariables(RatesTemplate):
        RNDLV = Float(-99.) # N loss rate leaves [kg ha-1 d-1]
        RNDST = Float(-99.) # N loss rate stems  [kg ha-1 d-1]
        RNDRT = Float(-99.) # N loss rate roots  [kg ha-1 d-1]
        
        RPDLV = Float(-99.) # P loss rate leaves [kg ha-1 d-1]
        RPDST = Float(-99.) # P loss rate stems  [kg ha-1 d-1]
        RPDRT = Float(-99.) # P loss rate roots  [kg ha-1 d-1]
        
        RKDLV = Float(-99.) # K loss rate leaves [kg ha-1 d-1]
        RKDST = Float(-99.) # K loss rate stems  [kg ha-1 d-1]
        RKDRT = Float(-99.) # K loss rate roots  [kg ha-1 d-1]
    
    class Parameters(ParamTemplate):
        RDRRTB = AfgenTrait() # rel. death rate of roots as a function of DVS [-; kg kg-1 d-1]
        RDRSTB = AfgenTrait() # rel. death rate of stems as a function of DVS [-; kg kg-1 d-1]
        
        RNFLV = Float(-99.) # residual N fraction in leaves [kg N kg-1 dry biomass]
        RNFST = Float(-99.) # residual N fraction in stems [kg N kg-1 dry biomass]
        RNFRT = Float(-99.) # residual N fraction in roots [kg N kg-1 dry biomass]
        RPFLV = Float(-99.) # residual P fraction in leaves [kg P kg-1 dry biomass]
        RPFST = Float(-99.) # residual P fraction in stems [kg P kg-1 dry biomass]
        RPFRT = Float(-99.) # residual P fraction in roots [kg P kg-1 dry biomass]
        RKFLV = Float(-99.) # residual K fraction in leaves [kg K kg-1 dry biomass]
        RKFST = Float(-99.) # residual K fraction in stems [kg K kg-1 dry biomass]
        RKFRT = Float(-99.) # residual K fraction in roots [kg K kg-1 dry biomass]
        
        
    def initialize(self, day, kiosk, cropdata):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PyWOFOST instance
        :param cropdata: dictionary with WOFOST cropdata key/value pairs
        :returns: the npk_losses_dying __call__()
        """

        self.params = self.Parameters(cropdata)
        self.rates   = self.RateVariables(kiosk,publish=["RNDLV","RNDST","RNDRT",
                                                         "RPDLV","RPDST","RPDRT",
                                                         "RKDLV","RKDST","RKDRT"])
        self.kiosk  = kiosk
        
    
    @prepare_rates
    def calc_rates(self, day):
        rates  = self.rates
        params = self.params
        
        DVS  = self.kiosk["DVS"]
        DRLV = self.kiosk["DRLV"] # death rate leaves [kg dry matter ha-1 d-1]
        
        WRT  = self.kiosk["WRT"]
        WST  = self.kiosk["WST"]
        
        
        DRRT = params.RDRRTB(DVS) * WRT
        DRST = params.RDRSTB(DVS) * WST
                
        rates.RNDLV = params.RNFLV * DRLV
        rates.RNDST = params.RNFST * DRST
        rates.RNDRT = params.RNFRT * DRRT
        
        rates.RPDLV = params.RPFLV * DRLV
        rates.RPDST = params.RPFST * DRST
        rates.RPDRT = params.RPFRT * DRRT

        rates.RKDLV = params.RKFLV * DRLV
        rates.RKDST = params.RKFST * DRST
        rates.RKDRT = params.RKFRT * DRRT
        

        
