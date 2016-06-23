# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from __future__ import print_function
import datetime
from copy import deepcopy

import numpy as np

from ..pydispatch import dispatcher
from ..traitlets import Float, Int, Instance, Enum, Unicode, Bool, AfgenTrait
from ..decorators import prepare_rates, prepare_states
from math import log10, sqrt, exp
from ..util import limit, Afgen, merge_dict
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, \
    SimulationObject
from .. import signals
from .. import exceptions as exc

def zeros(n):
    """Mimic np.zeros() by returning a list of zero floats of length n.
    """
    if isinstance(n, int):
        if n > 0:
            return [0.]*n
    
    msg = "zeros() should be called with positive integer, got: %s" % n
    raise ValueError(msg)

#-------------------------------------------------------------------------------
class WaterbalanceLayered(SimulationObject):
    """Waterbalance for freely draining soils under water-limited production.

    The purpose of the soil water balance calculations is to estimate the
    daily value of the soil moisture content. The soil moisture content
    influences soil moisture uptake and crop transpiration.
    
    The dynamic calculations are carried out in two sections, one for the 
    calculation of rates of change per timestep (= 1 day) and one for the
    calculation of summation variables and state variables. The water balance
    is driven by rainfall, possibly buffered as surface storage, and
    evapotranspiration. The processes considered are infiltration, soil water
    retention, percolation (here conceived as downward water flow from rooted
    zone to second layer), and the loss of water beyond the maximum root zone. 

    The textural profile of the soil is conceived as homogeneous. Initially the
    soil profile consists of two layers, the actually rooted  soil and the soil
    immediately below the rooted zone until the maximum rooting depth (soil and
    crop dependent). The extension of the root zone from initial rooting depth
    to maximum rooting depth is described in Root_Dynamics class. From the
    moment that the maximum rooting depth is reached the soil profile is
    described as a one layer system.
    
    The class WaterbalanceLayered is derived from WATFDGW.F90 in WOFOSTx.x
    (release March 2012)
    """
    # INTERNALS
    RDold  = Float(-99.)                # previous maximum rooting depth value
    RDMSLB = Float(-99.)                # max rooting depth soil layer boundary
    DSLR   = Float(-99.)                # Counter for Days-Dince-Last-Rain 
    RINold = Float(-99)                 # Infiltration rate of previous day
    
    XDEF   = Float(1000.0)              # maximum depth of groundwater (in cm)
    PFFC   = Float(2.0)                 # PF field capacity, Float(log10(200.))
    PFWP   = Float(log10(16000.))       # PF wilting point
    PFSAT  = Float(-1.0)                # PF saturation

    EquilTableLEN = Int(30)             # GW: WaterFromHeight, HeightFromAir
    MaxFlowIter = Int(50)
    
    # Fraction of non-infiltrating rainfall as function of storm size
    NINFTB = Instance(Afgen)

    # ------------------------------------------
    # ILaR: code taken from classic waterbalance
    # ------------------------------------------
    flag_crop_emerged = Bool(False)
    flag_crop_finished = Bool(False)
    # ------------------------------------------
    
    class Parameters(ParamTemplate):
        GW     = Int(-99)               # groundwater 
        ZTI    = Int(-99)               # inital depth of groundwater
        DD     = Int(-99)               # drainage depth
        RDMSOL = Float(-99.)            # not used for waterbalance

        # soil layers
        NSL      = Int(-99)             # from soildata, number
        SOIL_LAYERS= Instance(list)     # from soildata, layers (TSL+SOIL_GROUP_NO+typ)

        # crop parameters
        IAIRDU = Float(-99.)
        RDMCR  = Float(-99.)
        RDI    = Float(-99.)
        #RD     = Float(-99.)            # root_dynamics
        #RDM    = Float(-99.)            # root_dynamics
        # Site parameters
        IFUNRN = Float(-99.)            #Int(-99)?
        SSMAX  = Float(-99.)
        SSI    = Float(-99.)
        WAV    = Float(-99.)
        NOTINF = Float(-99.)
        SMLIM  = Float(-99.)
            
    class StateVariables(StatesTemplate):
        SM     = Float(-99.)
        SS     = Float(-99.)
        W      = Float(-99.)
        WI     = Float(-99.)
        WLOW   = Float(-99.)
        WBOT   = Float(-99.)
        WLOWI  = Float(-99.)
        WWLOW  = Float(-99.)            # W + WLOW
        WWLOWI = Float(-99.)
        WAVUPP = Float(-99.)
        WAVLOW = Float(-99.)
        RunOff = Float(-99.)
        #LOSS   = Float(-99.)
        #RIN    = Float(-99.)           # IN?
        
        # GW, groundwater
        WZ     = Float(-99.)
        WZI    = Float(-99.)
        ZT     = Float(-99.)
        WSUB0  = Float(-99.)
        WSUB   = Float(-99.)
        
        # layer
        ILR    = Int(-99)
        ILM    = Int(-99)

        # Summation variables 
        WTRAT  = Float(-99.)
        EVST   = Float(-99.)
        EVWT   = Float(-99.)
        TSR    = Float(-99.)
        RAINT  = Float(-99.)
        WDRT   = Float(-99.)
        TOTINF = Float(-99.)
        TOTIRR = Float(-99.)
        SUMSM  = Float(-99.)
        PERCT  = Float(-99.)
        LOSST  = Float(-99.)
        CRT    = Float(-99.)            # GW
        DRAINT = Float(0.)              # GW, always 0
        
        # Checksums for rootzone (RT) and total system (TT)
        WBALRT = Float(-99.)
        WBALTT = Float(-99.)

    class RateVariables(RatesTemplate):
        EVS   = Float(-99.)
        EVW   = Float(-99.)
        
        RIN   = Float(-99.)
        RIRR  = Float(0.)               # always 0
        PERC  = Float(-99.)
        LOSS  = Float(-99.)
        DW    = Float(-99.)
        DWLOW = Float(-99.)
        DWBOT = Float(-99.)             # GW
        DWSUB = Float(-99.)             # GW
        
        CR    = Float(-99.)
        DMAX  = Float(-99.)
        DZ    = Float(-99.)

        RAIN  = Float(-99.)             # weather
        WTRA  = Float(-99.)             # evapotranspiration
        WTRAL = Instance(np.ndarray)    # evapotranspiration

    class Layer(dict):                  # NOT USED as such ... in SOIL_LAYERS
        #soil_layers
        SOIL_GROUP_NO = Int(-99)
        TSL     = Float(-99.)           # thickness
    
        #soildata
        SMFCF   = Float(-99.)
        SM0     = Float(-99.)
        SMW     = Float(-99.)
        SOPE    = Float(-99.)           # soil surface conductivity
        KSUB    = Float(-99.)
        K0      = Float(-99.)
        #CRAIRC  = Float(-99.)          # not used for waterbalance
        CONTAB  = AfgenTrait()          # Conductivity from PF
        SMTAB   = AfgenTrait()          # Soil Moisture from PF
        
        #derived from soildata
        PFTAB   = AfgenTrait()          # Soil Moisture from PF
        MFPTAB  = AfgenTrait()          # Matrix Flux Potential
        WaterFromHeight = AfgenTrait()  # 
        HeightFromAir = AfgenTrait()    # 
        
        LBSL    = Float(-99.)           # lower boundary
        Wtop    = Float(-99.)           # layer weight factors
        Wpot    = Float(-99.)
        Wund    = Float(-99.)
    
        WC      = Float(-99.)
        WC0     = Float(-99.)
        WCW     = Float(-99.)
        WCFC    = Float(-99.)
        CondFC  = Float(-99.)
        CondK0  = Float(-99.)
        
        SM      = Float(-99.)
        DWC     = Float(-99.)
    
    def initialize(self, day, kiosk, cropdata, soildata, sitedata, missing):
        # Merge dictionaries in order to pass them to the Parameters class.
        # use merge_dict iso deepcopy?
        parvalues = deepcopy(cropdata)
        parvalues.update(soildata)
        parvalues.update(sitedata)
        parvalues.update(missing)
        parvalues["RDMSOL"] = parvalues["SOIL_LAYERS"][-1]["LBSL"]
        # Assign parameter values            
        self.params = self.Parameters(parvalues)
        
        p = self.params
        
        # ------ checks ------
        RD = p.RDI
        #classic: RDM = max(p.RDI, min(p.RDMSOL, p.RDMCR))
        RDM = max(p.RDI, p.RDMCR)
        
        #print "initialize WaterbalanceLayered NSL %i, GW %s, RD %f" % (p.NSL, p.GW, RD)
        
        if RD > RDM:
            msg = ("rooting depth %f exceeeds maximum rooting depth %f" %
                   (RD, RDM))
            raise exc.WaterBalanceError(msg)

        RDMFND = False
        for il in range (0, p.NSL):
            if abs(p.SOIL_LAYERS[il]['LBSL'] - RDM) < 0.01:
                RDMFND = True
                # layer boundary explicitly assigned to maximum rooting depth
                self.RDMSLB = p.SOIL_LAYERS[il]['LBSL']

        # also guarantees that RDM is within the layered part of the soil
        if not RDMFND:
            msg = ("Maximum rooting depth (RDM) %f does not coincide " +
                   "with a layer boundary in soil profile" % RDM)
            raise exc.WaterBalanceError(msg)
            
        # in case of groundwater the reference depth XDEF should be below the layered soil
        if p.GW:
            if self.XDEF <= p.SOIL_LAYERS[p.NSL-1]['LBSL']:
                msg = ("Reference depth XDEF (%f cm) must be below the " +
                       "bottom of the soil layers" % self.XDEF)
                raise exc.WaterBalanceError(msg)
        # --- end of checks ---

        # find deepest layer with roots
        ILR = p.NSL-1
        ILM = p.NSL-1
        for il in range(p.NSL-1, -1, -1):
            if (p.SOIL_LAYERS[il]['LBSL'] >= RD ):         ILR = il
            if (p.SOIL_LAYERS[il]['LBSL'] >= self.RDMSLB): ILM = il

        # calculate layer weight for RD-rooted layer and RDM-rooted layer
        self._layer_weights(RD, self.RDMSLB, ILR, ILM)
        # --- end of soil input section ---

        # save old rooting depth (for testing on growth in integration)
        self.RDold = RD

        SS = p.SSI                  # Initial surface storage

        # state variables set initially by self.StateVariables     
        W      = 0.0
        WAVUPP = 0.0
        WLOW   = 0.0
        WAVLOW = 0.0
        WBOT   = 0.0
        if p.GW:
            # calculate initial soil moisture
            ZT = limit(0.1, self.XDEF, p.ZTI)       # initial groundwater level
            if p.DD > 0.:                           # IDRAIN==1 ???
                ZT = max(ZT, p.DD)                  # corrected for drainage depth

            #for the soil layers   
            for il in range (0, p.NSL):
                HH = p.SOIL_LAYERS[il]['LBSL'] - p.SOIL_LAYERS[il]['TSL'] / 2.0   # depth at half-layer-height
                if p.SOIL_LAYERS[il]['LBSL'] - ZT < 0.0:
                    # layer is above groundwater ; get equilibrium amount from Half-Height pressure head
                    p.SOIL_LAYERS[il]['SM'] = p.SOIL_LAYERS[il]['SOILTYPE']['SMTAB'](log10(ZT-HH))
                elif p.SOIL_LAYERS[il]['LBSL']-p.SOIL_LAYERS[il]['TSL'] >= ZT:
                    # layer completely in groundwater
                    p.SOIL_LAYERS[il]['SM'] = p.SOIL_LAYERS[il]['SOILTYPE']['SM0']
                else:
                    # layer partly in groundwater
                    p.SOIL_LAYERS[il]['SM'] = (p.SOIL_LAYERS[il]['LBSL']-ZT)*p.SOIL_LAYERS[il]['SOILTYPE']['SM0'] \
                                                + p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](ZT-(p.SOIL_LAYERS[il]['LBSL']-p.SOIL_LAYERS[il]['TSL']) \
                                                ) / p.SOIL_LAYERS[il]['TSL']

            # calculate (available) water in rooted and potentially rooted zone
            # note that amounts WBOT below RDM (RDMSLB) are not available (below potential rooting depth)
            for il in range (0, p.NSL):
                W      += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wtop']
                WLOW   += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']
                WBOT   += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wund']
                # available water
                WAVUPP += (p.SOIL_LAYERS[il]['SM']-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) \
                            * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wtop']
                WAVLOW += (p.SOIL_LAYERS[il]['SM']-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) \
                            * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']

            # now various subsoil amonts
            WSUB0 = (self.XDEF-p.SOIL_LAYERS[p.NSL-1]['LBSL']) * p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['SM0']  # saturation
            if ZT > p.SOIL_LAYERS[p.NSL-1]['LBSL']:
                # groundwater below layered system
                WSUB = (self.XDEF-ZT)*p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['SM0'] \
                        + p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['WaterFromHeight'](ZT-p.SOIL_LAYERS[p.NSL-1]['LBSL'])
            else:
                # saturated subsoil
                WSUB = WSUB0

            # then amount of moisture below rooted zone
            WZ  = WLOW + WBOT + WSUB
            WZI = WZ

        else: # not GW
            # AVMAX - maximum available content of layer(s)
            # to get an even distribution of water in the rooted top if WAV is small.
            AVMAX  = np.zeros(p.NSL)
            TOPRED = Instance(Float)
            LOWRED = Instance(Float)
            SML    = Instance(Float)
            TOPLIM = 0.0
            LOWLIM = 0.0
            for il in range (0, ILM+1):
                # determine maximum content for this layer
                if il <= ILR:
                    # in the rooted zone a separate limit applies
                    SML = p.SMLIM
                    SML = limit(p.SOIL_LAYERS[il]['SOILTYPE']['SMW'], p.SOIL_LAYERS[il]['SOILTYPE']['SM0'], SML)
                    
                    # Check whether SMLIM is within boundaries
                    if p.IAIRDU == 1: # applicable only for flooded rice crops
                        SML = p.SOIL_LAYERS[il]['SOILTYPE']['SM0']
                        
                    # notify user of changes in SMLIM
                    if SML != p.SMLIM:
                        msg = "SMLIM not in valid range, changed from %f to %f."
                        self.logger.warn(msg % (p.SMLIM, SML))

                    AVMAX[il] = (SML-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) * p.SOIL_LAYERS[il]['TSL']   # available in cm  
                    # also if partly rooted, the total layer capacity counts in TOPLIM
                    # this means the water content of layer ILR is set as if it would be
                    # completely rooted. This water will become available after a little
                    # root growth and through numerical mixing each time step.
                    TOPLIM = TOPLIM + AVMAX[il]
                else:
                    # below the rooted zone the maximum is saturation (see code for WLOW in one-layer model)
                    # again the full layer capacity adds to LOWLIM. 
                    SML = p.SOIL_LAYERS[il]['SOILTYPE']['SM0']

                    AVMAX[il] = (SML-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) * p.SOIL_LAYERS[il]['TSL']   # available in cm  
                    LOWLIM += AVMAX[il]

            if p.WAV <= 0.0:
                # no available water
                TOPRED = 0.0
                LOWRED = 0.0
            elif p.WAV <= TOPLIM:
                # available water fits in layer(s) 1..ILR, these layers are rooted or almost rooted
                # reduce amounts with ratio WAV / TOPLIM
                TOPRED = p.WAV / TOPLIM
                LOWRED = 0.0
            elif p.WAV < TOPLIM+LOWLIM:
                # available water fits in potentially rooted layer
                # rooted zone is filled at capacity ; the rest reduced
                TOPRED = 1.0
                LOWRED = (p.WAV-TOPLIM) / LOWLIM
            else:
                # water does not fit ; all layers "full"
                TOPRED = 1.0
                LOWRED = 1.0

            # within rootzone
            for il in range (0, ILR+1):
                # Part of the water assigned to ILR may not actually be in the rooted zone, but it will
                # be available shortly through root growth (and through numerical mixing).
                p.SOIL_LAYERS[il]['SM'] = p.SOIL_LAYERS[il]['SOILTYPE']['SMW'] + AVMAX[il] * TOPRED /p.SOIL_LAYERS[il]['TSL']
                
                W      += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wtop']
                WLOW   += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']
                # available water
                WAVUPP += (p.SOIL_LAYERS[il]['SM']-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) \
                            * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wtop']
                WAVLOW += (p.SOIL_LAYERS[il]['SM']-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) \
                            * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']

            # between initial and maximum rooting depth. In case RDM is not a layer boundary (it should be!!)
            # layer ILM contains additional water in unrooted part. Only rooted part contributes to WAV.
            for il in range (ILR+1, ILM+1):
                p.SOIL_LAYERS[il]['SM'] = p.SOIL_LAYERS[il]['SOILTYPE']['SMW'] + \
                                        AVMAX[il] * LOWRED / p.SOIL_LAYERS[il]['TSL']
                
                WLOW += p.SOIL_LAYERS[il]['SM'] * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']
                # available water
                WAVLOW += (p.SOIL_LAYERS[il]['SM']-p.SOIL_LAYERS[il]['SOILTYPE']['SMW']) \
                            * p.SOIL_LAYERS[il]['TSL'] * p.SOIL_LAYERS[il]['Wpot']

            # below the maximum rooting depth
            for il in range (ILM+1, p.NSL):
                p.SOIL_LAYERS[il]['SM'] = p.SOIL_LAYERS[il]['SOILTYPE']['SMW']

            # set groundwater depth far away for clarity ; this prevents also
            # the root routine to stop root growth when they reach the groundwater
            ZT    = 999.0
            # init GW variables
            WZ    = 0
            WZI   = WZ
            WSUB0 = 0

        if p.GW:
            print("WATER LIMITED CROP PRODUCTION WITH GROUNDWATER")
        else:
            print("WATER LIMITED CROP PRODUCTION WITHOUT GROUNDWATER")
        print("=================================================")
        print(" fixed fraction       RDMso=%3.0f   NOTinf=%.3f" % (p.RDMSOL, p.NOTINF))
        print("                      SMLIM=%.3f   RDM=%3.0i  WAV=%3.0f  SSmax=%3.0f" % \
                                           (p.SMLIM, self.RDMSLB, p.WAV, p.SSMAX))
            
        # water content for each layer + a few fixed points often used
        for il in range (0, p.NSL):
            p.SOIL_LAYERS[il]['WC']     = p.SOIL_LAYERS[il]['SM']    * p.SOIL_LAYERS[il]['TSL']   # state variable
            p.SOIL_LAYERS[il]['WC0']    = p.SOIL_LAYERS[il]['SOILTYPE']['SM0']   * p.SOIL_LAYERS[il]['TSL']
            p.SOIL_LAYERS[il]['WCW']    = p.SOIL_LAYERS[il]['SOILTYPE']['SMW']   * p.SOIL_LAYERS[il]['TSL']
            p.SOIL_LAYERS[il]['WCFC']   = p.SOIL_LAYERS[il]['SOILTYPE']['SMFCF'] * p.SOIL_LAYERS[il]['TSL']
            p.SOIL_LAYERS[il]['CondFC'] = 10.0**p.SOIL_LAYERS[il]['SOILTYPE']['CONTAB'](self.PFFC)
            p.SOIL_LAYERS[il]['CondK0'] = 10.0**p.SOIL_LAYERS[il]['SOILTYPE']['CONTAB'](self.PFSAT)
            #print "layer %i: TSL %f WC %f SM %f SMfcf %f SMsat %f" % (il, \
            # p.SOIL_LAYERS[il]['TSL'], p.SOIL_LAYERS[il]['WC'], p.SOIL_LAYERS[il]['SM'], \
            # p.SOIL_LAYERS[il]['SOILTYPE']['SMFCF'], p.SOIL_LAYERS[il]['SOILTYPE']['SM0'])
            print("layer  %i  %3.1f cm: SM0=%.3f SMFC=%.3f SMW=%.3f" % (il, p.SOIL_LAYERS[il]['TSL'], \
                p.SOIL_LAYERS[il]['SOILTYPE']['SM0'], p.SOIL_LAYERS[il]['SOILTYPE']['SMFCF'], \
                p.SOIL_LAYERS[il]['SOILTYPE']['SMW']))

        # rootzone and subsoil water
        WI    = W
        WLOWI = WLOW
        WWLOW = W + WLOW
        SM = W / RD         # p.SOIL_LAYERS[p.NSL-1]['LBSL']

        # soil evaporation, days since last rain @ToDo: <= or <?
        self.DSLR = 1.0
        if p.SOIL_LAYERS[0]['SM'] <= (p.SOIL_LAYERS[0]['SOILTYPE']['SMW'] + \
                                      0.5*(p.SOIL_LAYERS[0]['SOILTYPE']['SMFCF'] - \
                                           p.SOIL_LAYERS[0]['SOILTYPE']['SMW'])):
            self.DSLR=5.0
        self.RINold = 0.    # RIN is used in calc_rates before it is set, so keep the RIN as RINold?
        self.NINFTB = Afgen([0.0,0.0, 0.5,0.0, 1.5,1.0, 0.0,0.0, 0.0,0.0, \
                             0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0])
        
        # INITIAL STATES
        self.states = self.StateVariables(kiosk, publish=["SM","W","WLOW","RAINT","WTRAT"],
            SM=SM, SS=SS, ZT=ZT, WZ=WZ, WZI=WZI, WSUB0=WSUB0,
            W=W, WI=WI, WLOW=WLOW, WLOWI=WLOWI, WWLOW=WWLOW, WWLOWI=WWLOW,
            WAVUPP=WAVUPP, WAVLOW=WAVLOW, ILR=ILR, ILM=ILM,
            WTRAT=0.0, EVST=0.0, EVWT=0.0, TSR=0.0,
            RAINT=0.0, WDRT=0.0, TOTINF=0.0, TOTIRR=0.0, SUMSM=0.0,
            PERCT=0.0, LOSST=0.0, RunOff=0.0, IN=0.0,
            CRT=0.0, DRAINT=0.0, WBOT=0.0, WSUB=0.0,
            WBALRT=-999.0, WBALTT=-999.0)
        
        self.rates = self.RateVariables(kiosk)
        
        self.kiosk = kiosk

        # ------------------------------------------
        # ILaR: code taken from classic waterbalance
        # ------------------------------------------
        # Connect to CROP_EMERGED/CROP_FINISH signals for water balance to
        # search for crop transpiration values
        self._connect_signal(self._on_CROP_EMERGED, signals.crop_emerged)
        self._connect_signal(self._on_CROP_FINISH, signals.crop_finish)
        # ------------------------------------------
    #---------------------------------------------------------------------------

    @prepare_rates
    def calc_rates(self, day, drv): # layered
        s = self.states
        p = self.params
        r = self.rates

        DELT = 1
        RD = self._determine_rooting_depth()
        if RD != self.RDold:
            msg = "Rooting depth changed unexpectedly"
            raise RuntimeError(msg)

        #print "calc_rates WaterbalanceLayered NSL %i, GW %s, RD %f" % (p.NSL, p.GW, RD)
        
        # conductivities and Matric Flux Potentials for all layers
        PF            = np.zeros(p.NSL)
        Conductivity  = np.zeros(p.NSL)
        MatricFluxPot = np.zeros(p.NSL)
        EquilWater    = np.zeros(p.NSL)
        for il in range (0, p.NSL):
            PF[il] = p.SOIL_LAYERS[il]['SOILTYPE']['PFTAB'](p.SOIL_LAYERS[il]['SM'])     
            #print "layer %i: PF %f SM %f" % (il, PF[il], p.SOIL_LAYERS[il]['SM'])
            Conductivity[il]  = 10.0**p.SOIL_LAYERS[il]['SOILTYPE']['CONTAB'](PF[il])
            MatricFluxPot[il] = p.SOIL_LAYERS[il]['SOILTYPE']['MFPTAB'](PF[il])
            
            if p.GW:  # equilibrium amounts
                if p.SOIL_LAYERS[il]['LBSL'] < s.ZT:
                    # groundwater below layer
                    EquilWater[il] = p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](s.ZT-p.SOIL_LAYERS[il]['LBSL']+\
                                                                                      p.SOIL_LAYERS[il]['TSL']) \
                                   - p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](s.ZT-p.SOIL_LAYERS[il]['LBSL'])
                elif p.SOIL_LAYERS[il]['LBSL']-p.SOIL_LAYERS[il]['TSL'] < s.ZT:
                    # groundwater in layer
                    EquilWater[il] = p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](s.ZT-p.SOIL_LAYERS[il]['LBSL']+\
                                                                                      p.SOIL_LAYERS[il]['TSL']) \
                                   + (p.SOIL_LAYERS[il]['LBSL']-s.ZT) * p.SOIL_LAYERS[il]['SOILTYPE']['SM0']
                else: # groundwater above layer
                    EquilWater[il] = p.SOIL_LAYERS[il]['WC0']

        # ------------------------------------------
        # ILaR: code taken from classic waterbalance
        # ------------------------------------------
        # Rainfall rate, not from kiosk???
        r.RAIN = drv.RAIN

        # Transpiration and maximum soil and surfacewater evaporation rates
        # are calculated by the crop Evapotranspiration module. 
        # However, if the crop is not yet emerged then set TRA=0 and use
        # the potential soil/water evaporation rates directly because there is
        # no shading by the canopy.
        if self.flag_crop_emerged is True:
            r.WTRA  = self.kiosk["TRA"]
            r.WTRAL = self.kiosk["TRALY"]
            EVWMX   = self.kiosk["EVWMX"]
            EVSMX   = self.kiosk["EVSMX"]
        else:
            r.WTRA  = 0.
            r.WTRAL = np.zeros(20)
            EVWMX   = drv.E0
            EVSMX   = drv.ES0
        # ------------------------------------------

        # actual evaporation rates ...
        r.EVW = 0.
        r.EVS = 0.
        # ... from surface water if surface storage more than 1 cm, ...
        if s.SS > 1.:
            r.EVW = EVWMX
        else:    # ... else from soil surface
            if self.RINold >= 1.:       # RIN not set, must be RIN from previous 'call'
                r.EVS  = EVSMX
                self.DSLR = 1.
            else:
                self.DSLR += 1.
                EVSMXT = EVSMX*(sqrt(self.DSLR) - sqrt(self.DSLR-1.))
                r.EVS  = min(EVSMX, EVSMXT+self.RINold)

        # preliminary infiltration rate
        if s.SS <= 0.1:   # without surface storage
            if p.IFUNRN==0.: RINPRE = (1.-p.NOTINF)*r.RAIN + r.RIRR + s.SS/DELT
            if p.IFUNRN==1.: RINPRE = (1.-p.NOTINF*self.NINFTB(r.RAIN))*r.RAIN + r.RIRR + s.SS/DELT
        else:
            # with surface storage, infiltration limited by SOPE (topsoil)
            AVAIL  = s.SS + (r.RAIN * (1.-p.NOTINF) + r.RIRR - r.EVW) * DELT
            RINPRE = min(p.SOIL_LAYERS[0]['SOILTYPE']['SOPE']*DELT, AVAIL) / DELT

        # maximum flow at Top Boundary of each layer
        # ------------------------------------------
        # DOWNWARD flows are calculated in two ways, 
        # (1) a "dry flow" from the matric flux potentials
        # (2) a "wet flow" under the current layer conductivities and downward gravity. 
        # Clearly, only the dry flow may be negative (=upward). The dry flow accounts for the large 
        # gradient in potential under dry conditions (but neglects gravity). The wet flow takes into 
        # account gravity only and will dominate under wet conditions. The maximum of the dry and wet 
        # flow is taken as the downward flow, which is then further limited in order the prevent 
        # (a) oversaturation and (b) water content to decrease below field capacity.
        #
        # UPWARD flow is just the dry flow when it is negative. In this case the flow is limited
        # to a certain fraction of what is required to get the layers at equal potential, taking
        # into account, however, the contribution of an upward flow from further down. Hence, in 
        # case of upward flow from the groundwater, this upward flow in propagated upward if the 
        # suction gradient is sufficiently large.

        EVflow  = np.zeros(p.NSL+1)         # 1 more
        FlowMX  = np.zeros(p.NSL+1)         # 1 more
        Flow    = np.zeros(p.NSL+1)         # 1 more
        LIMWET  = np.zeros(p.NSL)
        LIMDRY  = np.zeros(p.NSL)
        
        for il in range (0, p.NSL):
            p.SOIL_LAYERS[il]['DWC'] = 0.0  # water change

        # first get flow through lower boundary of bottom layer
        if p.GW:
            # the old capillairy rise routine is used to estimate flow to/from the groundwater
            # note that this routine returns a positive value for capillairy rise and a negative
            # value for downward flow, which is the reverse from the convention in WATFDGW.
            if s.ZT >= p.SOIL_LAYERS[p.NSL-1]['LBSL']:
                # groundwater below the layered system ; call the old capillairty rise routine
                # the layer PF is allocated at 1/3 * TSL above the lower boundary ; this leeds
                # to a reasonable result for groundwater approaching the bottom layer 
                SubFlow = self._SUBSOL(PF[p.NSL-1], \
                                       s.ZT-p.SOIL_LAYERS[p.NSL-1]['LBSL']+p.SOIL_LAYERS[p.NSL-1]['TSL']/3.0, \
                                       p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['CONTAB'])

                if SubFlow >= 0.0:
                    # capillairy rise is limited by the amount required to reach equilibrium:
                    # step 1. calculate equilibrium ZT for all air between ZT and top of layer
                    EqAir   = s.WSUB0 - s.WSUB + (p.SOIL_LAYERS[p.NSL-1]['WC0'] - p.SOIL_LAYERS[p.NSL-1]['WC'])
                    # step 2. the grouindwater level belonging to this amount of air in equilibrium
                    ZTeq1   = (p.SOIL_LAYERS[p.NSL-1]['LBSL'] - p.SOIL_LAYERS[p.NSL-1]['TSL']) \
                            + p.SOIL_LAYERS[il]['SOILTYPE']['HeightFromAir'](EqAir)
                    # step 3. this level should normally lie below the current level
                    #         (otherwise there should not be capillairy rise)
                    #         in rare cases however, due to the use of a mid-layer height
                    #         in subroutine SUBSOL, a deviation could occur
                    ZTeq2   = max(s.ZT, ZTeq1)
                    # step 4. calculate for this ZTeq2 the equilibrium amount of water in the layer
                    WCequil = p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](ZTeq2-p.SOIL_LAYERS[p.NSL-1]['LBSL']+\
                                                                               p.SOIL_LAYERS[p.NSL-1]['TSL']) \
                            - p.SOIL_LAYERS[il]['SOILTYPE']['WaterFromHeight'](ZTeq2-p.SOIL_LAYERS[p.NSL-1]['LBSL'])
                    # step5. use this equilibrium amount to limit the upward flow
                    FlowMX[p.NSL] = -1.0 * min(SubFlow, max(WCequil-p.SOIL_LAYERS[p.NSL-1]['WC'], 0.0)/DELT)
                else:
                    # downward flow ; air-filled pore space of subsoil limits downward flow
                    AirSub = (s.ZT-p.SOIL_LAYERS[p.NSL-1]['LBSL'])*p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['SM0'] \
                           - p.SOIL_LAYERS[p.NSL-1]['SOILTYPE']['WaterFromHeight'](s.ZT-p.SOIL_LAYERS[p.NSL-1]['LBSL'])
                    FlowMX[p.NSL] = min(abs(SubFlow), max(AirSub, 0.0)/DELT)
            else:
                # groundwater is in the layered system ; no further downward flow
                FlowMX[p.NSL] = 0.0
        else:   # not GW
            # Bottom layer conductivity limits the flow. Below field capacity there is no 
            # downward flow, so downward flow through lower boundary can be guessed as
            FlowMX[p.NSL] = max(p.SOIL_LAYERS[p.NSL-1]['CondFC'], Conductivity[p.NSL-1])

        # drainage
        r.DMAX = 0.0

        for il in range (p.NSL-1, -1, -1):
            # if this layers contains maximum rootig depth and if rice, downward water loss is limited
            if p.IAIRDU==1 and il==s.ILM:
                FlowMX[il+1] = 0.05 * p.SOIL_LAYERS[il]['CondK0']

            # limiting DOWNWARD flow rate 
            # == wet conditions: the soil conductivity is larger
            #    the soil conductivity is the flow rate for gravity only
            #    this limit is DOWNWARD only
            # == dry conditions: the MFP gradient
            #    the MFP gradient is larger for dry conditions
            #    allows SOME upward flow
            if il==0:
                LIMWET[il] = p.SOIL_LAYERS[0]['SOILTYPE']['SOPE']
                LIMDRY[il] = 0.0
            else:
                # same soil type
                if p.SOIL_LAYERS[il-1]['SOILTYPE']['SOIL_GROUP_NO'] == p.SOIL_LAYERS[il]['SOILTYPE']['SOIL_GROUP_NO']:
                    # flow rate estimate from gradient in Matric Flux Potential
                    LIMDRY[il] = 2.0 * (MatricFluxPot[il-1]-MatricFluxPot[il])/(p.SOIL_LAYERS[il-1]['TSL'] + p.SOIL_LAYERS[il]['TSL'])
                    if LIMDRY[il] < 0.0:
                        # upward flow rate ; amount required for equal water content is required below
                        MeanSM = (p.SOIL_LAYERS[il-1]['WC']+p.SOIL_LAYERS[il]['WC']) / (p.SOIL_LAYERS[il-1]['TSL'] + p.SOIL_LAYERS[il]['TSL'])
                        EqualPotAmount = p.SOIL_LAYERS[il-1]['WC'] - p.SOIL_LAYERS[il-1]['TSL'] * MeanSM  # should be negative like the flow
                        
                else:   # different soil types
                    # iterative search to PF at layer boundary (by bisection)
                    PF1  = PF[il-1]
                    PF2  = PF[il]
                    MFP1 = MatricFluxPot[il-1]
                    MFP2 = MatricFluxPot[il]
                    for i in range (0, self.MaxFlowIter):
                        PFx = (PF1 + PF2) / 2.0
                        Flow1 = 2.0 * (+ MFP1 - p.SOIL_LAYERS[il-1]['MFPTAB'](PFx)) / p.SOIL_LAYERS[il-1]['TSL']
                        Flow2 = 2.0 * (- MFP2 + p.SOIL_LAYERS[ il ]['MFPTAB'](PFx)) / p.SOIL_LAYERS[il]['TSL']
                        if abs(Flow1-Flow2) < TinyFlow: # sufficient accuracy
                            break
                        elif abs(Flow1) > abs(Flow2):
                            # flow in layer 1 is larger ; PFx must shift in the direction of PF1
                            PF2 = PFx
                        elif abs(Flow1) < abs(Flow2):
                            # flow in layer 2 is larger ; PFx must shift in the direction of PF2
                            PF1 = PFx

                    if i >= self.MaxFlowIter:
                        msg = "LIMDRY flow iteration failed"
                        raise RuntimeError(msg)
                    
                    LIMDRY[il] = (Flow1 + Flow2) / 2.0
                    if LIMDRY[il] < 0.0:
                        # upward flow rate ; amount required for equal potential is required below
                        Eq1 = -p.SOIL_LAYERS[il]['WC']
                        Eq2 = 0.0
                        for i in range (0, self.MaxFlowIter):
                            EqualPotAmount = (Eq1 + Eq2) / 2.0
                            SM1 = (p.SOIL_LAYERS[il-1]['WC'] - EqualPotAmount) / p.SOIL_LAYERS[il-1]['TSL']
                            SM2 = (p.SOIL_LAYERS[ il ]['WC'] + EqualPotAmount) / p.SOIL_LAYERS[il]['TSL']
                            PF1 = p.SOIL_LAYERS[il-1]['SMTAB'](SM1)
                            PF2 = p.SOIL_LAYERS[ il ]['SMTAB'](SM2)
                            if abs(Eq1-Eq2) < TinyFlow:     # sufficient accuracy
                                break
                            elif PF1 > PF2: # suction in top layer larger; absolute amount should be larger
                                Eq2 = EqualPotAmount
                            else:           # suction in bottom layer larger; absolute amount should be reduced
                                Eq1 = EqualPotAmount

                        if i >= self.MaxFlowIter:
                            msg = "Limiting amount iteration failed"
                            raise RuntimeError(msg)

                # the limit under wet conditions in a unit gradient
                LIMWET[il] = (p.SOIL_LAYERS[il-1]['TSL'] + p.SOIL_LAYERS[il]['TSL']) \
                            / (p.SOIL_LAYERS[il-1]['TSL']/Conductivity[il-1] + p.SOIL_LAYERS[il]['TSL']/Conductivity[il])

            FlowDown = True
            if LIMDRY[il] < 0.0:
                # upward flow (negative !) is limited by fraction of amount required for equilibrium
                FlowMax = max(LIMDRY[il], EqualPotAmount * UpwardFlowLimit)
                if il > 0:
                    # upward flow is limited by amount required to bring target layer at equilibrium/field capacity
                    if p.GW:  # soil does not drain below equilibrium with groundwater
                        FCequil = max(p.SOIL_LAYERS[il-1]['WCFC'], EquilWater[il-1])
                    else:   # free drainage
                        FCequil = p.SOIL_LAYERS[il-1]['WCFC']

                    TargetLimit = r.WTRAL[il-1] + (FCequil-p.SOIL_LAYERS[il-1]['WC'])/DELT
                    if TargetLimit > 0.0:
                        # target layer is "dry": below field capacity ; limit upward flow
                        FlowMax = max(FlowMax, -1.0*TargetLimit)
                        # there is no saturation prevention since upward flow leads to a decrease ofp.SOIL_LAYERS[il]['WC']
                        # instead flow is limited in order to prevent a negative water content
                        FlowMX[il] = max(FlowMax, FlowMX[il+1] + r.WTRAL[il] - p.SOIL_LAYERS[il]['WC']/DELT)
                        FlowDown = False
                    elif p.GW:
                        # target layer is "wet": above field capacity, since gravity is
                        # neglected in the matrix potential model, upward flow tends to be
                        # overestyimated in wet conditions. With groundwater the profile
                        # can get filled with water from above and upward flow is set to zero here.
                        FlowMX[il] = 0.0
                        FlowDown   = False
                    else:
                        # target layer is "wet": above field capacity, no groundwater
                        # free drainage model implies that upward flow is rejected here
                        # instead, downward flow is enabled. This guarantees that, if all
                        # layers are above field capacity, the free drainage model applies.
                        FlowDown = True

            if FlowDown:
                # maximum downward flow rate (LIMWET is always a positive number)
                FlowMax = max(LIMDRY[il], LIMWET[il])
                # this prevents saturation of layer il
                # maximum top boundary flow is bottom boundary flow plus saturation deficit plus sink
                FlowMX[il] = min(FlowMax, FlowMX[il+1] + (p.SOIL_LAYERS[il]['WC0']-p.SOIL_LAYERS[il]['WC'])/DELT + r.WTRAL[il])

        # adjustment of infiltration rate to prevent saturation
        r.RIN = min(RINPRE, FlowMX[0])

        # contribution of layers to soil evaporation in case of drought upward flow is allowed
        EVSL    = np.zeros(p.NSL)
        EVSL[0] = min(r.EVS, (p.SOIL_LAYERS[0]['WC'] - p.SOIL_LAYERS[0]['WCW'])/DELT + r.RIN - r.WTRAL[0])
        
        EVrest  = r.EVS - EVSL[0]
        for il in range (1, p.NSL):
            Available = max(0.0, (p.SOIL_LAYERS[il]['WC'] - p.SOIL_LAYERS[il]['WCW'])/DELT - r.WTRAL[il])
            if Available >= EVrest:
                EVSL[il] = EVrest
                EVrest   = 0.0
                break
            else:
                EVSL[il] = Available
                EVrest   = EVrest - Available

        # reduce evaporation if entire profile becomes airdry
        # there is no evaporative flow through lower boundary of layer NSL
        r.EVS -= EVrest

        #! evaporative flow (taken positive) at layer boundaries
        EVflow[0] = r.EVS
        for il in range (1, p.NSL):
            EVflow[il] = EVflow[il-1] - EVSL[il-1]
        EVflow[p.NSL] = 0.0 

        # limit downward flows not to get below field capacity / equilibrium content         
        Flow[0] = r.RIN - EVflow[0]
        #print "Flow %i: %f = RIN %f - EVflow %f, sm: %f" % \
        # (0, Flow[0], r.RIN, EVflow[0], p.SOIL_LAYERS[0]['SM'])
        
        for il in range (0, p.NSL):
            if p.GW:  # soil does not drain below equilibrium with groundwater
                WaterLeft = max(p.SOIL_LAYERS[il]['WCFC'], EquilWater[il])
            else:    # free drainage
                WaterLeft = p.SOIL_LAYERS[il]['WCFC']

            MXLOSS     = (p.SOIL_LAYERS[il]['WC']-WaterLeft)/DELT   # maximum loss 
            #print "   MXLOSS %i: %f = WC %f - WaterLeft %f, sm: %f" % \
            # (il, MXLOSS, p.SOIL_LAYERS[il]['WC'], WaterLeft, p.SOIL_LAYERS[il]['SM'])
            Excess     = max(0.0, MXLOSS + Flow[il] - r.WTRAL[il])  # excess of water (positive)
            #print "   Excess %i: %f = max(0.0, MXLOSS %f + Flow %f - WTRAL %f)" % \
            # (il, Excess, MXLOSS, Flow[il], r.WTRAL[il])
            Flow[il+1] = min(FlowMX[il+1], Excess - EVflow[il+1])   # negative (upward) flow is not affected
            #print "Flow   %i: %f = min(FlowMX %f, Excess %f - EVflow %f)" %
            # (il+1, Flow[il+1], FlowMX[il+1], Excess, EVflow[il+1])
        
            # rate of change
            p.SOIL_LAYERS[il]['DWC'] = Flow[il] - Flow[il+1] - r.WTRAL[il]
            #print "layer %i: DWC %f=%f-%f-%f, sm: %f" % (il,\
            #    p.SOIL_LAYERS[il]['DWC'], Flow[il], Flow[il+1], r.WTRAL[il], \
            #    p.SOIL_LAYERS[il]['SM'])

        # Percolation and Loss.
        # Equations were derived from the requirement that in the same layer, above and below 
        # depth RD (or RDM), the water content is uniform. Note that transpiration above depth 
        # RD (or RDM) may require a negative percolation (or loss) for keeping the layer uniform. 
        # This is in fact a numerical dispersion. After reaching RDM, this negative (LOSS) can be 
        # prevented by choosing a layer boundary at RDM.
        if s.ILR < s.ILM:
            # layer ILR is devided into rooted part (where the sink is) and a below-roots part
            # The flow in between is PERC
            f1 = p.SOIL_LAYERS[s.ILR]['Wtop']   # 1-f1 = Wpot
            r.PERC = (1.0-f1) * (Flow[s.ILR]-r.WTRAL[s.ILR]) + f1 * Flow[s.ILR+1]

            # layer ILM is divided as well ; the flow in between is LOSS
            f2 = p.SOIL_LAYERS[s.ILM]['Wpot']
            f3 = 1.0 - f2		        # f3 = Wund
            r.LOSS = f3 * Flow[s.ILM] + f2 * Flow[s.ILM+1]
        elif s.ILR == s.ILM:
            # depths RD and RDM in the same soil layer: there are three "sublayers":
            # - the rooted sublayer with fraction f1
            # - between RD and RDM with fraction f2
            # - below RDM with fraction f3
            # PERC goes from 1->2, LOSS from 2->3
            # PERC and LOSS are calculated in such a way that the three sublayers have equal SM
            f1 = p.SOIL_LAYERS[s.ILR]['Wtop']
            r.PERC = (1.0-f1) * (Flow[s.ILR]-r.WTRAL[s.ILR]) + f1 * Flow[s.ILR+1]
            
            f2 = p.SOIL_LAYERS[s.ILM]['Wpot']
            f3 = 1.0 - f1 - f2
            r.LOSS = f3 * (Flow[s.ILM]-r.WTRAL[s.ILM]) + (1.0-f3) * Flow[s.ILM+1]
        else:
            msg = "Internal_1"
            raise RuntimeError(msg)

        # rates of change in amounts of moisture W and WLOWI
        r.DW    = - r.WTRAL.sum() - r.EVS - r.PERC + r.RIN
        r.DWLOW = r.PERC - r.LOSS
        #print "DW    %f= -WTRAL %f -EVS %f -PERC %f +RIN %f" % (r.DW, r.WTRAL.sum(), r.EVS, r.PERC, r.RIN)
        #print "DWLOW %f= PERC %f - LOSS %f" % (r.DWLOW, r.PERC, r.LOSS)

        if p.GW:  # groundwater influence
            r.DWBOT = r.LOSS - Flow[p.NSL]
            r.DWSUB = Flow[p.NSL]
            #print "DWBOT %f= LOSS %f - Flow %f" % (r.DWBOT, r.LOSS, Flow[p.NSL])
            #print "DWSUB %f= Flow %f" % (r.DWSUB, Flow[p.NSL])
        #msg = '\n'.join(['%s = %s' % (k,v) for k,v in self.kiosk.iteritems()]) 
        #print(msg)
                
    #---------------------------------------------------------------------------

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        p = self.params
        r = self.rates
        
        #print "integrate WaterbalanceLayered NSL %i, GW %s" % (p.NSL, p.GW)
        
        DELT = 1    # zou weg kunnen en dan ook overal DELT opruimen
        
        #!-----------------------------------------------------------------------
        #! integrals of the water balance:  summation and state variables
        #!-----------------------------------------------------------------------
        #! amount of water in soil layers ; soil moisture content
        for il in range (0, p.NSL):
            p.SOIL_LAYERS[il]['WC'] += p.SOIL_LAYERS[il]['DWC'] * DELT
            p.SOIL_LAYERS[il]['SM']  = p.SOIL_LAYERS[il]['WC'] / p.SOIL_LAYERS[il]['TSL']
            #print "layer %i: WC +%f=%f sm: %f" % (il, p.SOIL_LAYERS[il]['DWC'], p.SOIL_LAYERS[il]['WC'], p.SOIL_LAYERS[il]['SM'])

        # totals
        s.WTRAT  += r.WTRAL.sum()*DELT      # transpiration
        s.EVWT   += r.EVW*DELT              # evaporation from surface water layer and/or soil
        s.EVST   += r.EVS*DELT
        s.RAINT  += r.RAIN*DELT             # rainfall, irrigation and infiltration
        s.TOTINF += r.RIN*DELT
        s.TOTIRR += r.RIRR*DELT
        
        SSPRE  = s.SS + (r.RAIN+r.RIRR-r.EVW-r.RIN)*DELT
        s.SS   = min(SSPRE, p.SSMAX)        # surface storage
        s.TSR += (SSPRE-s.SS)               # surface runoff
        
        # amounts of water
        s.W      += r.DW*DELT               # in rooted zone
        # classic test if negative: if W < 0.0: EVST += W ; w = 0.0
        s.WAVUPP += r.DW*DELT
        s.WLOW   += r.DWLOW*DELT            # in unrooted, lower part of rootable zone
        s.WAVLOW += r.DWLOW*DELT
        
        s.WWLOW = s.W + s.WLOW              # total in the whole rootable zone
        s.WBOT += r.DWBOT * DELT            # and in layered soil below RDM
        
        # percolation from rootzone ; interpretation depends on mode
        if p.GW:                    # flow is either percolation or capillairy rise
            if r.PERC > 0.0:
               s.PERCT += r.PERC*DELT
            else:
               s.CRT   -= r.PERC*DELT
        else:                       # flow is always called percolation
            s.PERCT += r.PERC*DELT
            s.CRT    = 0.0

        # loss of water by flow from the potential rootzone
        s.LOSST += r.LOSS*DELT

        #----------------------------------------------
        # change of rootzone subsystem boundary
        #----------------------------------------------
        # calculation of amount of soil moisture in new rootzone
        RD = self._determine_rooting_depth()                # ???
        if (RD-self.RDold) > 0.001:
            # roots have grown find new values ; overwrite W, WLOW, WAVUPP, WAVLOW, WBOT
            s.ILR = p.NSL-1
            s.ILM = p.NSL-1
            for il in range (p.NSL-1, -1, -1):
                if (p.SOIL_LAYERS[il]['LBSL'] >= RD ):         s.ILR = il
                if (p.SOIL_LAYERS[il]['LBSL'] >= self.RDMSLB): s.ILM = il

            self._layer_weights(RD, self.RDMSLB, s.ILR, s.ILM)

            WOLD = s.W
            s.W    = 0.0
            s.WLOW = 0.0
            s.WBOT = 0.0
            s.WAVUPP = 0.0
            s.WAVLOW = 0.0

            # get W and WLOW and available water amounts
            for il in range (0, p.NSL):
                s.W    += p.SOIL_LAYERS[il]['WC'] * p.SOIL_LAYERS[il]['Wtop']
                s.WLOW += p.SOIL_LAYERS[il]['WC'] * p.SOIL_LAYERS[il]['Wpot']
                s.WBOT += p.SOIL_LAYERS[il]['WC'] * p.SOIL_LAYERS[il]['Wund']

                s.WAVUPP += (p.SOIL_LAYERS[il]['WC']-p.SOIL_LAYERS[il]['WCW']) * p.SOIL_LAYERS[il]['Wtop']
                s.WAVLOW += (p.SOIL_LAYERS[il]['WC']-p.SOIL_LAYERS[il]['WCW']) * p.SOIL_LAYERS[il]['Wpot']
                
            WDR     = s.W - WOLD    # water added to root zone by root growth, in cm
            s.WDRT += WDR           # total water addition to rootzone by root growth
            
            self.RDold = RD         # save RD for which layer contents have been determined

        s.SM = s.W/RD               # mean soil moisture content in rooted zone
        s.SUMSM += s.SM*DELT        # calculating mean soil moisture content over growing period

        #----------------------------------------------
        # groundwater level
        #----------------------------------------------
        if p.GW:  # with groundwater influence
            s.WSUB += r.DWSUB * DELT            # subsoil between soil layers and reference plane
            s.WZ    = s.WLOW + s.WBOT + s.WSUB  # amount of water below rooted zone
            
            # find groundwater level
            ZTfound = False
            for il in range (p.NSL-1, -1, -1):
                if il==p.NSL-1:
                    AirSub = s.WSUB0 - s.WSUB
                    if AirSub > 0.01:
                        # groundwater is in subsoil which is not completely saturated
                        s.ZT = min(p.SOIL_LAYERS[p.NSL-1]['LBSL'] + p.SOIL_LAYERS[il]['SOILTYPE']['HeightFromAir'](AirSub), self.XDEF)
                        ZTfound = True
                        break

                if p.SOIL_LAYERS[il]['SM'] < 0.999 * p.SOIL_LAYERS[il]['SOILTYPE']['SM0']:
                    # groundwater is in this layer
                    s.ZT = p.SOIL_LAYERS[il]['LBSL'] - p.SOIL_LAYERS[il]['TSL'] \
                        + min(p.SOIL_LAYERS[il]['TSL'], p.SOIL_LAYERS[il]['SOILTYPE']['HeightFromAir'](p.SOIL_LAYERS[il]['WC0']-p.SOIL_LAYERS[il]['WC']))
                    ZTfound = True
                    break

            if not ZTfound: # entire system saturated
                s.ZT = 0.0

        # quick-and-dirty: do once some printing on leaves died
        if self.kiosk.get("WLV",-9) == 0 and s.WBALRT == -999:
            self.finalize(day)
    #---------------------------------------------------------------------------

    @prepare_states
    def finalize(self, day):
        s = self.states
        p = self.params

        # Checksums waterbalance for rootzone (WBALRT) and whole system (WBALTT)
        # ---
        # GW: WZI/WZ iso WLOWI/WLOW, + CRT, DRAINT iso LOSST, LOSST part of WZ?
        if p.GW:
            s.WBALRT = s.TOTINF + s.CRT + s.WI + s.WDRT \
                       - s.EVST - s.WTRAT - s.PERCT - s.W
            s.WBALTT = p.SSI + s.RAINT + s.TOTIRR + s.WI + s.WZI \
            - s.W - s.WZ - s.WTRAT - s.EVWT - s.EVST - s.TSR - s.DRAINT - s.SS
        else:
            # mean water content rooting zone during crop growth and total
            # water content of the potentially rooted zone at end of simulation
            #MWC = SUMSM/max(1.,REAL (MOD((365+IDHALT-IDEM),365)))
            #TWE = W+WLOW

            s.WBALRT = s.TOTINF + s.WI + s.WDRT \
                       - s.EVST - s.WTRAT - s.PERCT - s.W
            s.WBALTT = p.SSI + s.RAINT + s.TOTIRR + s.WI - s.W + s.WLOWI - \
                       s.WLOW - s.WTRAT - s.EVWT - s.EVST - s.TSR - s.LOSST - s.SS

            print("\n       WATER BALANCE WHOLE SYSTEM (1 DIMENS. COLUMN ; cm)")
            print(" init max root zone   %5.1f  final max root zone  %5.1f    change: %5.1f" % (s.WWLOWI, s.WWLOW, s.WWLOW-s.WWLOWI))
            print("  init surf storage   %5.1f  final surf storage   %5.1f    change: %5.1f" % (p.SSI, s.SS, s.SS-p.SSI))
            print("         irrigation   %5.1f  evap water surface   %5.1f"                  % (s.TOTIRR, s.EVWT))
            print("           rainfall   %5.1f  evap soil surface    %5.1f"                  % (s.RAINT, s.EVST))
            print("                             transpiration        %5.1f  to atmos: %5.1f" % (s.WTRAT, s.EVWT+s.EVST+s.WTRAT))
            print("                             surface runoff       %5.1f"                  % (s.TSR))
            print("                             lost to deep soil    %5.1f"                  % (s.LOSST))
            print("  TOTAL INIT + IN     %5.1f  TOTAL FINAL + OUT    %5.1f  checksum: %5.1f" % (s.WWLOWI+p.SSI+s.TOTIRR+s.RAINT, \
                                                                                                s.WWLOW+s.SS+s.EVWT+s.EVST+s.WTRAT+s.TSR+s.LOSST, \
                                                                                                s.WBALTT))
            
            print("\n               WATER BALANCE ROOT ZONE")
            print(" initial water stock  %5.1f  final water stock    %5.1f"                  % (s.WI, s.W))
            print("        infiltration  %5.1f  evap soil surface    %5.1f"                  % (s.TOTINF, s.EVST))
            print(" added by root growth %5.1f      transpiration    %5.1f"                  % (s.WDRT, s.WTRAT))
            print("                                   percolation    %5.1f"                  % (s.PERCT))
            print(" TOTAL INIT + IN      %5.1f  FINAL + OUT          %5.1f  checksum: %5.1f" % (s.WI+s.TOTINF+s.WDRT, \
                                                                                                s.W+s.EVST+s.WTRAT+s.PERCT, \
                                                                                                s.WBALRT))
            print("\n")
            
        if abs(s.WBALRT) > 0.0001 or abs(s.WBALTT) > 0.0001:
            msg = "Error in layered waterbalance!"
            raise RuntimeError(msg)

        # Run finalize on the subSimulationObjects
        SimulationObject.finalize(self, day)

    def _determine_rooting_depth(self):
        """Determines appropriate use of the rooting depth (RD)
        """
        p = self.params
        
        if self.flag_crop_emerged is False and \
           self.flag_crop_finished is False:
            
            # before crop starts hold RD at initial value RDI
            RD = p.RDI
            
        elif self.flag_crop_emerged is True and \
             self.flag_crop_finished is False:
            
            # In cropping season
            RD = self.kiosk["RD"]

        elif self.flag_crop_emerged is True and \
             self.flag_crop_finished is True:
            
            # Crop finished
            if  "RD" in self.kiosk:
                # Only happens at the final simulation cycle when value for
                # SM still has to be computed
                RD = self.kiosk["RD"]
            else:
                RD = p.RDI
        else:
            msg = "Unable to determine rooting depth in Waterbalance Layered"
            raise RuntimeError(msg)
        
        return RD

    def _layer_weights(self, RD, RDM, ILR, ILM):
        """Calculate weight factors for rooted- and sub-layer calculations
        """
        # RD    the rooting depth
        # ILR   deepest layer containing roots
        # ILM   deepest layer that will contain roots
        # RDM   max rooting depth aligned to soil layer boundary
        
        # NSL   number of layers
        # TSL   the layerthickness
        # LBSL  the Lower Boundaries of the NSL soil layers
        # Wtop  weights for contribution to rootzone
        # Wpot  weights for contribution to potentially rooted zone
        # Wund  weights for contribution to never rooted layers

        p   = self.params
        s   = self.states
        NSL = p.NSL
        sl  = p.SOIL_LAYERS
        
        #print "---\nlayer_weights NSL: %i ILR: %i ILM: %i RD: %f RDM: %f\n---" % \
        # (NSL, ILR, ILM, RD, RDM)

        for il in range (0, NSL):
            # rooted layer
            if (il < ILR):                  
                sl[il]['Wtop'] = 1.0
                sl[il]['Wpot'] = 0.0
                sl[il]['Wund'] = 0.0
                
            # partly rooted
            elif (il == ILR and il < ILM):  # at the end fully rooted
                sl[il]['Wtop'] = 1.0 - (sl[il]['LBSL']-RD) / sl[il]['TSL']
                sl[il]['Wpot'] = 1.0 - sl[il]['Wtop']
                sl[il]['Wund'] = 0.0
            elif (il == ILR and il == ILM): # at the end partly rooted
                sl[il]['Wtop'] = 1.0 - (sl[il]['LBSL']-RD) / sl[il]['TSL']
                sl[il]['Wund'] = (sl[il]['LBSL']-RDM) / sl[il]['TSL']
                sl[il]['Wpot'] = 1.0 - sl[il]['Wund'] - sl[il]['Wtop']
                
            # not rooted
            elif (il < ILM):                # at the end fully rooted
                sl[il]['Wtop'] = 0.0
                sl[il]['Wpot'] = 1.0
                sl[il]['Wund'] = 0.0
            elif (il == ILM):               # at the end partly rooted
                sl[il]['Wtop'] = 0.0
                sl[il]['Wund'] = (sl[il]['LBSL']-RDM) / sl[il]['TSL']
                sl[il]['Wpot'] = 1.0 - sl[il]['Wund']
                
            # never rooted
            else:                           
                sl[il]['Wtop'] = 0.0
                sl[il]['Wpot'] = 0.0
                sl[il]['Wund'] = 1.0
                
            #msg = 'layer %i: %f %f weights top %f pot %f und %f' % \
            #       (il, sl[il]['TSL'], sl[il]['LBSL'], \
            #        sl[il]['Wtop'], sl[il]['Wpot'], sl[il]['Wund'])
            #print msg

    def _SUBSOL(self, PF, D, CONTAB):
        """SUBSOL...
        """
        
        DEL = np.zeros(4)
        PFGAU = np.zeros(12)
        HULP = np.zeros(12)
        CONDUC = np.zeros(12)
        
        ELOG10 = 2.302585
        LOGST4 = 2.518514
        
        START  = np.array([0, 45, 170, 330])
        PFSTAN = np.array([0.705143, 1.352183, 1.601282, 1.771497, 2.031409, \
                           2.192880, 2.274233, 2.397940, 2.494110])
        PGAU = np.array([0.1127016654, 0.5, 0.8872983346])
        WGAU = np.array([0.2777778, 0.4444444, 0.2777778])
        
        # calculation of matric head and check on small pF
        PF1 = PF
        D1  = D
        MH  = exp(ELOG10*PF1)
        
        if PF1 <= 0.:
            # in case of small matric head
            K0   = exp(ELOG10 * CONTAB(-1))
            FLOW = K0 * (MH/D-1)
        else:
            IINT = 0
            
            # number and width of integration intervals
            for I1 in range (0, 4):
                if I1 <= 2:
                    DEL[I1] = min(START[I1+1],MH)-START[I1]
                elif I1 == 3:
                    DEL[I1] = PF1-LOGST4
                    
                if DEL[I1] <= 0:
                    break;
                
                IINT += 1
    
            # preparation of three-point Gaussian integration
            for I1 in range (0, IINT):
                for I2 in range (0, 3):
                    I3 = (3*I1) + I2
                    
                    if I1 == IINT-1:
                        # the three points in the last interval are calculated
                        if IINT <= 3:
                            PFGAU[I3] = log10 (START[IINT-1] + PGAU[I2] * DEL[IINT-1])
                        elif IINT == 4:
                            PFGAU[I3] = LOGST4 + PGAU[I2] * DEL[IINT-1]
                    else:
                        PFGAU[I3] = PFSTAN[I3]

                    # variables needed in the loop below
                    CONDUC[I3] = exp(ELOG10 * CONTAB(PFGAU[I3]))
                    HULP[I3]   = DEL[I1] * WGAU[I2] * CONDUC[I3]
                    
                    if I3 > 8:
                        HULP[I3] = HULP[I3] * ELOG10 * exp(ELOG10 * PFGAU[I3]);
    
            # 15.5 setting upper and lower limit
            FU = 1.27
            FL = -1 * exp(ELOG10 * CONTAB(PF1))
            if MH <= D1: FU = 0
            if MH >= D1: FL = 0
            if MH != D1:
                # Iteration loop
                IMAX = 3 * IINT;
                
                for I1 in range (0, 15):
                    FLW = (FU+FL)/2
                    DF  = (FU-FL)/2
                    if DF < 0.01 and DF/abs(FLW) < 0.1: 
                        break

                    Z = 0
                    for I2 in range (0, IMAX):
                        Z += HULP[I2]/(CONDUC[I2] + FLW)

                    if Z >= D1: FL = FLW
                    if Z <= D1: FU = FLW
    
            FLOW = (FU+FL)/2
                    
        return FLOW

    def _on_CROP_EMERGED(self):
        self.flag_crop_emerged = True
        
    def _on_CROP_FINISH(self):
        self.flag_crop_finished = True
