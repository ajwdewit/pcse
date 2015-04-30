'''
LINTUL2_N with adaptation to N-Limitation
'''
from pcse.base_classes import SimulationObject, ParamTemplate, StatesTemplate,\
    RatesTemplate
from pcse.traitlets import Float, AfgenTrait
from pcse.decorators import prepare_rates, prepare_states
from pcse.lintul import lintul3lib
from pcse.lintul.lintul3lib import notNull, INSW, REAAND
from numpy.ma.core import exp
from start import Lintul3Model


class Lintul3(SimulationObject):
    """
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

************************************************************************
*   LINTUL3 is an extended version of LINTUL1 (the version of LINTUL   *
*          for optimal growth conditions) and LINTUL2 includes a simple*
*          water balance for studying effects of drought. LINTUL2-N    *
*          includes N-limitation on crop growth. The latter program is *
*          called LINTUL3.                                             *
*          Test version for spring wheat, using parameters for spring  *
*          wheat for Flevoland                                         *
************************************************************************
"""
    # Parameters, rates and states which are relevant at the main crop
    # simulation level
    class Parameters(ParamTemplate):
        DOYEM = Float(-99.)
        DRATE = Float(-99.)
        DVSDR = Float(-99.)
        DVSNLT= Float(-99.)
        DVSNT = Float(-99.)
        FNTRT = Float(-99.)
        FRNX  = Float(-99.)
        IRRIGF= Float(-99.)
        K     = Float(-99.)
        LAICR = Float(-99.)
        LRNR  = Float(-99.)
        LSNR  = Float(-99.)
        LUE   = Float(-99.)
        NFRLVI= Float(-99.)
        NFRRTI= Float(-99.)
        NFRSTI= Float(-99.)
        NLAI  = Float(-99.)
        NLUE  = Float(-99.)
        NMAXSO= Float(-99.)
        NPART = Float(-99.)
        NSLA  = Float(-99.)
        RDRSHM= Float(-99.)
        RGRL  = Float(-99.)
        RNFLV = Float(-99.)
        RNFRT = Float(-99.)
        RNFST = Float(-99.)
        ROOTDM= Float(-99.)
        RRDMAX= Float(-99.)
        SLAC  = Float(-99.)
        TBASE = Float(-99.)
        TCNT  = Float(-99.)
        TRANCO= Float(-99.)
        TSUMAG= Float(-99.)
        TSUMAN= Float(-99.)
        TSUMMT= Float(-99.)
        WCAD  = Float(-99.)
        WCFC  = Float(-99.)
        WCI   = Float(-99.)
        WCST  = Float(-99.)
        WCSUBS= Float(-99.)
        WCWET = Float(-99.)
        WCWP  = Float(-99.)
        WMFAC = Float(-99.)

        FERTAB = AfgenTrait()
        FLVTB  = AfgenTrait()
        FRTTB  = AfgenTrait()
        FSOTB  = AfgenTrait()
        FSTTB  = AfgenTrait()
        NMXLV  = AfgenTrait()
        NRFTAB = AfgenTrait()
        PHOTTB = AfgenTrait()
        RDRT   = AfgenTrait()
        SLACF  = AfgenTrait()
        


        
    class StateVariables(StatesTemplate):
        TSUM  = Float(-99.)
        LAI   = Float(-99.)
        ANLV  = Float(-99.)
        ANST  = Float(-99.)
        ANRT  = Float(-99.)
        ANSO  = Float(-99.)
        NUPTT = Float(-99.)
        TNSOIL= Float(-99.)
        NLOSSL= Float(-99.)
        NLOSSR= Float(-99.)
        WLVG  = Float(-99.)
        WLVD  = Float(-99.)
        WST   = Float(-99.)
        WSO   = Float(-99.)
        WRT   = Float(-99.)
        ROOTD = Float(-99.)
        GTSUM = Float(-99.)
        WDRT  = Float(-99.)
        CUMPAR= Float(-99.)
        WA    = Float(-99.)
        TEXPLO= Float(-99.)
        TEVAP = Float(-99.)
        TTRAN = Float(-99.)
        TRUNOF= Float(-99.)
        TIRRIG= Float(-99.)
        TRAIN = Float(-99.)
        TDRAIN= Float(-99.)
#         
#                 
#         def __getattr__(self, name)
#             tryState = name[]
#             if (name.starts_with('r') and hasattr(self, name)):
#                 
#             return StatesTemplate.__getattribute__(self, *args, **kwargs)

        @classmethod
        def listIntegratedStates(cls):
            return sorted([a for a in cls.__dict__ if isinstance(getattr(cls, a), Float) and not a.startswith('_')])





    class InitialValues(object):
        
        def __init__(self, parameters):
            # Read initial states
            self.ROOTDI= 0.1
            self.TSUMI = 0.0
            self.WLVGI = 2.4
            self.WRTLI = 3.6
            self.WSOI  = 0.0
            self.WSTI  = 0.0    
            
            # initial calculations
            DVSI= self.TSUMI / parameters.TSUMAN
            SLACFI = parameters.SLACF(DVSI)
            ISLA   = parameters.SLAC * SLACFI
            
            # Initial amount of water present in the rooted depth at the start of
            # the calculations, based on the initial water content (in mm).
            self.WAI  = 1000. * self.ROOTDI * parameters.WCI
            
            # Initial amount of N (g/m2) in leaves, stem, roots, and storage organs.
            self.ANLVI = parameters.NFRLVI * self.WLVGI
            self.ANSTI = parameters.NFRSTI * self.WSTI
            self.ANRTI = parameters.NFRRTI * self.WRTLI
            self.ANSOI = 0.
            
            #   Initial LAI.
            self.LAII   = self.WLVGI * ISLA
            
        
        
    @classmethod
    def defineRateVariables(cls):
        '''
        ensure a rate for each state
        :param cls: outer class of RateVariables  and StateVariables
        '''
        attributes = {}
        for s in cls.StateVariables.listIntegratedStates():
            attributes['r' + s] = Float(-99.)

        return type('RateVariables', (RatesTemplate, ), attributes)
        
        
        
    initialValues = None
    ratesClass = None
    DSLR  = 0 
    FSHMOD = 0.0       
    OUTPUT_VARS = sorted( ["TIME", "WAI", "DVS", "TSUM", "TAGBM", "WST", "WLVG", "WLVD", "WSO", "LAI", "NTAC", "WRT", 
               "GTSUM", "CBALAN", "TRANRF", "NNI", "SLA", "FRACT", "FRTWET", "FLVT", "FSTT", "FSOT", 
               "RWLVG", "RWST", "RWRT", "RWSO", "CUMPAR", "LUECAL", "NUPTT", "TTRAN", "TEVAP", "PEVAP", 
               "NBALAN", "WATBAL", "NUPTR", "TNSOIL", "NDEMTO", "RNSOIL", "FERTN", "FERTNS", "WA", 
               "TIRRIG", "TRAIN", "TEXPLO", "TRUNOF", "TDRAIN"])
 
    
    def __init__(self, day, kiosk, *args, **kwargs):
        self.find_subroutines()
        
        super(Lintul3, self).__init__(day, kiosk, *args, **kwargs)
        
        
    
    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.kiosk  = kiosk
        self.params = self.Parameters(parvalues)

        # Read initial states
        init                = self.InitialValues(self.params)
        self.initialValues  = init
        initialStates       = dict((a, 0.0) for a in self.StateVariables.listIntegratedStates())
        
        
        # Initialize state variables
        initialStates["TSUM"]   = init.TSUMI
        initialStates["LAI"]    = init.LAII
        initialStates["ANLV"]   = init.ANLVI
        initialStates["ANST"]   = init.ANSTI
        initialStates["ANRT"]   = init.ANRTI
        initialStates["ANSO"]   = init.ANSOI
        initialStates["WLVG"]   = init.WLVGI
        initialStates["WST"]    = init.WSTI
        initialStates["WSO"]    = init.WSOI
        initialStates["WRT"]    = init.WRTLI
        initialStates["ROOTD"]  = init.ROOTDI
        initialStates["WA"]     = init.WAI

        self.states = self.StateVariables(kiosk, publish=[], **initialStates)
                
#         self.rates  = self.RateVariables(kiosk)
        self.ratesClass = self.defineRateVariables()
        self.rates  = self.ratesClass(kiosk)

        
                
    def find_subroutines(self):
        self.astro  = lintul3lib.ASTRO
        self.deathl = lintul3lib.DEATHL
        self.drunir = lintul3lib.DRUNIR
        self.evaptr = lintul3lib.EVAPTR
        self.gla    = lintul3lib.GLA
        self.growth = lintul3lib.GROWTH
        self.ndemnd = lintul3lib.NDEMND
        self.nnindx = lintul3lib.NNINDX
        self.noptm  = lintul3lib.NOPTM
        self.ntrans = lintul3lib.NTRANS
        self.ntrloc = lintul3lib.NTRLOC
        self.penman = lintul3lib.PENMAN
        self.relgr  = lintul3lib.RELGR
        self.rnld   = lintul3lib.RNLD
        self.rnusub = lintul3lib.RNUSUB
        self.subdvs = lintul3lib.SUBDVS
        self.subpar = lintul3lib.SUBPAR
        
            
    
    
    @prepare_rates
    def calc_rates(self, day, drv):
        # dynamic calculations
        p = self.params
        r = self.rates
        s = self.states
        i = self.initialValues
        
        TIME =  day.timetuple().tm_yday
        DELT = 1 # ???
        
        __out__ = (TIME==180)
        
        # Variables supplied by the weather system
        LAT              = drv.LAT
        RDD              = drv.IRRAD # ???
        TMMN             = drv.TMIN
        TMMX             = drv.TMAX
        VP               = drv.VAP  / 10 # hPa --> kPa correction from Cabo Weather
        WN               = drv.WIND
        RAIN             = drv.RAIN * 10 # cm  --> mm CORRECTION FOR NON-STANDARD cm in CABO-WEATHER
                
        # ----------Emergence, Temperature sum and Developmental stages----------*
        
        DTR    = RDD/1.E+6  # Actual daily total global radiation (DTR, J m-2 d-1, the factor 1.E06 converts J into MJ)
        DAVTMP = 0.5 * (TMMN + TMMX)
##      TTSUM  = p.TSUMAN + p.TSUMMT     ## obsolete
        
            
        # **********************Calling Subroutines******************************
        #  Calling the subroutine for calculating the astrological daylength.
        DAYL = self.astro(TIME, LAT)
        
        # Calling the subroutine for converting TSUM to the developmental stage.
        DVS = self.subdvs(TIME, p.DOYEM, s.TSUM, p.TSUMAN, p.TSUMMT)
        
        # Calling the subroutine for translocatable N in leaves, stem, roots and
        # storage organs.
        ATNLV, ATNST, ATNRT, ATN = self.ntrloc(s.ANLV, s.ANST, s.ANRT, s.WLVG, s.WST, s.WRT, 
                                               p.RNFLV, p.RNFST, p.RNFRT, p.FNTRT)
        
        # *Total vegetative biomass.
        TBGMR = s.WLVG + s.WST
        

        # N concentration (g N g-1 DM) of the leaves, stem, roots and storage
        # organs
##      NFLV = s.ANLV/ notNull(s.WLVG)  ## obsolete
##      NFST = s.ANST/ notNull(s.WST)   ## obsolete
##      NFRT = s.ANRT/ notNull(s.WRT)   ## obsolete
##      NFSO = s.ANSO/ notNull(s.WSO)   ## obsolete
        
        # Total N in green matter of the plant.
        NUPGMR = s.ANLV + s.ANST
        
        # ---------------Fertilizer application---------------------------------*
        FERTN  = p.FERTAB(TIME)
        NRF    = p.NRFTAB(TIME)
        

        
        #  Relative death rate of leaves due to N stress.
        RDRNS  = 0.03
        
        # Total leaf weight.
        WLV    = s.WLVG + s.WLVD
        
        #  Water content in the rootzone
        WC  = 0.001* s.WA /notNull(s.ROOTD)
        
        #  Relative death rate of roots.
        RDRRT = 0.03
        
        
        # ----------------------------------------------------------------------*
        
##      NTAG   = s.ANLV + s.ANST + s.ANSO   ## obsolete
        
        # Relative death rate of leaves due to senescence/ageing.
        RDRTMP = p.RDRT(DAVTMP)
        
        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        
        NMAXLV = p.NMXLV(DVS)
        
        # Photoperiodic effect.
        PHOTPF  = INSW (s.TSUM-p.TSUMAN, p.PHOTTB(DAYL), 1.)
        
        # * Total above ground biomass
        TAGBM = WLV + s.WST + s.WSO ## --> output
        
        # N supply to the storage organs.
        NSUPSO = INSW (DVS - p.DVSNT, 0., ATN / p.TCNT)
        EMERG  = REAAND(TIME - p.DOYEM + 1., WC - p.WCWP)* INSW(-s.LAI, 1., 0.)
        
        # Calling the subroutine for Potential evaporation and transpiration.
        # RLWN, NRADC, PENMRC, PENMD, ...
        PEVAP, PTRAN = self.penman(DAVTMP, VP, DTR, s.LAI, WN)
        NFGMR  = NUPGMR / notNull(TBGMR)
        
        # *Average residual N concentration.
        NRMR   = (s.WLVG * p.RNFLV + s.WST * p.RNFST) / notNull(TBGMR)
        
        #  Nitrogen uptake limiting factor at low moisture conditions in the
        #  rooted soil layer before anthesis. After anthesis there is no
        #  uptake from the soil anymore.
        NLIMIT = INSW(DVS - p.DVSNLT, INSW(WC - p.WCWP, 0., 1.), 0.0)
        FERTNS = FERTN * NRF
        
        # -------- Growth rates and dry matter production of plant organs-------*
        #  Biomass partitioning functions under (water and nitrogen)non-stressed
        #  situations
        FRTWET = p.FRTTB( DVS )
        FLVT   = p.FLVTB( DVS )
        FSTT   = p.FSTTB( DVS)
        FSOT   = p.FSOTB( DVS )
        
        # * For calculation of LUE (it can be removed once after checking, however, 
        #   I have put these new additions in these two (c and N) balances sections
        PAR    = DTR * 0.50
        DTEFF  = max ( 0., DAVTMP - p.TBASE )
        RTSUM  = DTEFF * EMERG
        
        # Calling the subroutine for actual rates of evaporation and
        # transpiration.
#         WCCR, EVAP, TRAN, FR, AVAILF, self.DSLR = self.evaptr(PEVAP, PTRAN, s.ROOTD, s.WA, 
        EVAP, TRAN, self.DSLR = self.evaptr(PEVAP, PTRAN, s.ROOTD, s.WA, 
                                            p.WCAD, p.WCWP, p.WCFC, p.WCWET, p.WCST, 
                                            p.TRANCO, DELT, p.WMFAC, RAIN, self.DSLR)
        NMAXST = p.LSNR * NMAXLV
        NMAXRT = p.LRNR * NMAXLV
        
        #  Soil N supply (g N m-2 d-1) through mineralization.
        RTMIN  = 0.10 * EMERG * NLIMIT
        RROOTD = min(p.RRDMAX * INSW(WC - p.WCWP, 0., 1.) * EMERG,  p.ROOTDM - s.ROOTD)
#         NTAC   = NTAG/notNull(TAGBM)
#         RTSH   = s.WRT / TAGBM
        RTSUMP = RTSUM * PHOTPF
        
        # Calling the subroutine for rates of drainage, runoff and irrigation.
        DRAIN, RUNOFF, IRRIG = self.drunir(RAIN, EVAP, TRAN, p.IRRIGF, p.DRATE, 
                                           DELT, s.WA, s.ROOTD, p.WCFC, p.WCST, p.WMFAC)
        
        #  Rainfall interception by the crop (strange use; left out)
        #  RNINTC = min( RAIN, 0.25*LAI )
        
        #  Exploration of water in soil when roots grow downward.
        EXPLOR = 1000. * RROOTD * p.WCSUBS
        
        # *   Growth reduction function for water stress(actual trans/potential)
        TRANRF = TRAN / notNull(PTRAN)
        
        # Calling the subroutine for maximum nitrogen concentration of leaves
        # and stem.
        NOPTLV, NOPTST = self.noptm(p.FRNX, NMAXLV, NMAXST)
        
        # Maximum N content in the plant.
        NOPTS = NOPTST * s.WST
        NOPTL = NOPTLV * s.WLVG
        RWA = (RAIN+EXPLOR+IRRIG)-(RUNOFF+TRAN+EVAP+DRAIN)
        NOPTMR = (NOPTL+ NOPTS)/notNull(TBGMR)
        
        # Calling the subroutine for Nitrogen Nutrition Index.
        NNI = self.nnindx(TIME, p.DOYEM, EMERG, NFGMR, NRMR, NOPTMR)
        
        # Calling the subroutine for relative modification for root and shoot
        # allocation.
        self.FSHMOD, FRT, FLV, FST, FSO = self.subpar(p.NPART, TRANRF, NNI, FRTWET, FLVT, FSTT, FSOT, self.FSHMOD)
        
        # Calling the subroutine for total growth rate.
        # PARINT, ... 
        GTOTAL = self.growth(TIME, p.DOYEM, DTR, p.K, p.NLUE, s.LAI, p.LUE, TRANRF, NNI)
        
        #  Water-Nitrogen stress factor
##      RNW = min(TRANRF, NNI)  ## obsolete
        
        # **-------------Functions and parameters for rice----------------------*
        
        #     Specific Leaf area(m2/g).
        SLA = p.SLAC * p.SLACF(DVS) * exp(-p.NSLA * (1.-NNI))
##      FRACT  = FLV+ FST+ FSO + FRT        ## obsolete
##      LUECAL = GTOTAL / PAR               ## obsolete
        
        # Calling the subroutine for relative death rate of leaves.
#         RDRDV, RDRSH, RDR, DLV, DLVS, DLVNS, DLAIS, DLAINS, DLAI = self.deathl(TIME, p.DOYEM, s.TSUM, p.TSUMAG, RDRTMP, 
        DLV, DLAI = self.deathl(TIME, p.DOYEM, s.TSUM, p.TSUMAG, RDRTMP, 
                                p.RDRSHM, s.LAI, p.LAICR, s.WLVG, RDRNS, NNI, SLA)
        
        # ** Leaf growth and LAI.
        GLV    = FLV * GTOTAL
        
        # Calling the subroutine for daily increase of leaf area index.
        GLAI = self.gla(TIME, p.DOYEM, DTEFF, i.LAII, p.RGRL, DELT, SLA, s.LAI, GLV, p.NLAI, WC, p.WCWP, DVS, TRANRF, NNI)
        
        # Calling the subroutine for N loss due to death of leaves and roots.
        DRRT, RNLDLV, RNLDRT = self.rnld(DVS, s.WRT, RDRRT, p.RNFLV, DLV, p.RNFRT, p.DVSDR)
        
        # Net rate of change of Leaf area.
        RLAI   = GLAI - DLAI
        
        # Calling the subroutine for relative growth rate of roots, leaves, stem
        # and storage organs.
        RWLVG, RWRT, RWST, RWSO = self.relgr(TIME, p.DOYEM, EMERG, 
                                             i.WLVGI, i.WRTLI, i.WSTI, i.WSOI, 
                                             GTOTAL, FLV, FRT, FST, FSO, DLV, DRRT, DELT)
        
        # Calling the subroutine for N demand of leaves, roots and stem storage
        # organs.
        NDEML, NDEMS, NDEMR, NDEMSO = self.ndemnd(NMAXLV, NMAXST, NMAXRT, p.NMAXSO, 
                                                  s.WLVG, s.WST, s.WRT, s.WSO, RWLVG, RWST, RWRT, RWSO, 
                                                  s.ANLV, s.ANST, s.ANRT, s.ANSO, p.TCNT, DELT)
        
        # ****************SOIL NITROGEN SUPPLY***********************************
        #     Total Nitrogen demand.
        NDEMTO = max(0.0, (NDEML + NDEMS + NDEMR))
        
        # Rate of N uptake in grains.
        RNSO =  min(NDEMSO, NSUPSO)
        
        # --------------------Fertilizer N---------------------------------------
        #     Total nutrient uptake.
        NUPTR = (max(0., min(NDEMTO, s.TNSOIL))* NLIMIT ) / DELT * INSW(TIME - p.DOYEM, 0., 1.)
        
        # Calling the subroutine for N translocated from leaves, stem, and roots.
        RNTLV, RNTST, RNTRT = self.ntrans(RNSO, ATNLV, ATNST, ATNRT, ATN)
        
        # Calling the subroutine to compute the partitioning of the total
        # N uptake rate (NUPTR) over the leaves, stem and roots.
        RNULV, RNUST, RNURT = self.rnusub(TIME, p.DOYEM, EMERG, NDEML, NDEMS, NDEMR, NUPTR, NDEMTO, 
                                          i.ANLVI, i.ANSTI, i.ANRTI, DELT)
        
        #  Change in inorganic N in soil as function of fertilizer
        #  input, soil N mineralization and crop uptake.
        RNSOIL = FERTNS/DELT -NUPTR + RTMIN
        RNST = RNUST-RNTST
        RNRT = RNURT-RNTRT-RNLDRT
        
        # *Rate of change of N in organs
        RNLV = RNULV-RNTLV-RNLDLV


        # *  Carbon, Nitrogen, Water balance
        CBALAN = (s.GTSUM + (i.WRTLI + i.WLVGI + i.WSTI + i.WSOI)     # @UnusedVariable
                         - (WLV + s.WST + s.WSO + s.WRT + s.WDRT))   

        NBALAN = (s.NUPTT + (i.ANLVI + i.ANSTI + i.ANRTI + i.ANSOI)   # @UnusedVariable
                         - (s.ANLV + s.ANST + s.ANRT + s.ANSO + s.NLOSSL + s.NLOSSR)) 

        WATBAL = (s.WA + (s.TRUNOF + s.TTRAN + s.TEVAP + s.TDRAIN)    # @UnusedVariable
                         - (i.WAI + s.TRAIN + s.TEXPLO +s.TIRRIG)) 
                 

            
        r.rTSUM   = RTSUMP
        r.rLAI    = RLAI  
        r.rANLV   = RNLV  
        r.rANST   = RNST  
        r.rANRT   = RNRT  
        r.rANSO   = RNSO  
        r.rNUPTT  = NUPTR 
        r.rTNSOIL = RNSOIL
        r.rNLOSSL = RNLDLV
        r.rNLOSSR = RNLDRT
        r.rWLVG   = RWLVG 
        r.rWLVD   = DLV   
        r.rWST    = RWST  
        r.rWSO    = RWSO  
        r.rWRT    = RWRT  
        r.rROOTD  = RROOTD
        r.rGTSUM  = GTOTAL
        r.rWDRT   = DRRT  
        r.rCUMPAR = PAR   
        r.rWA     = RWA   
        r.rTEXPLO = EXPLOR
        r.rTEVAP  = EVAP  
        r.rTTRAN  = TRAN  
        r.rTRUNOF = RUNOFF
        r.rTIRRIG = IRRIG 
        r.rTRAIN  = RAIN  
        r.rTDRAIN = DRAIN
        
        
        if (TIME == 1):
            for v in self.OUTPUT_VARS:
                try:
                    if locals().has_key(v):
                        value = locals()[v]
                    else:
                        value = getattr(self.states, v)
                except:
                    self.OUTPUT_VARS.remove(v)
            
            print "day\t",
            for v in self.OUTPUT_VARS:
                print "%s\t" % v,
            print
                        
        print "%s\t" % day,
        for v in self.OUTPUT_VARS:
            try:
                    if locals().has_key(v):
                        value = locals()[v]
                    else:
                        value = getattr(self.states, v)
            except:
                value =-99999
            print "%f\t" % (value),
            
        print
        
        
# finish conditions
#         if (KEEP == 1):
#             #    ---------------------------------------------------------------------*
#             #    Run control.
#             if (DVS >  2.01) or (TSUM > TTSUM): 
#                 TERMNL = true


    
    
    @prepare_states
    def integrate(self, day):
        rates = self.rates
        states = self.states
        delta = 1.
        
        TIME =  day.timetuple().tm_yday
        __out__ = (TIME==180)
        
        for s in states.listIntegratedStates():
            rate = getattr(rates, 'r' + s)
            state = getattr(states, s)
            newvalue = state + delta * rate
            setattr(states, s, newvalue)
#             if __out__:
#                 print "%s: %f = %f + %f" % (s, newvalue, state, rate)
            

        # crop stage before integration
#      crop_stage = self.pheno.get_variable("STAGE")

        # Phenology
#      self.pheno.integrate(day)
    


    @prepare_states
    def finalize(self, day):
       
        SimulationObject.finalize(self, day)
        
        
        
    @prepare_rates
    def updateRates(self):
        self.rates.rTSUM = 1
        
if (__name__ == "__main__"):

    sim = Lintul3Model.start(1987)
    l = sim.crop
    sim.run(182)

