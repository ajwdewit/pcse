'''
LINTUL2_N with adaptation to N-Limitation
'''
from pcse.base_classes import SimulationObject, ParamTemplate
from pcse.traitlets import Float, AfgenTrait, Instance
from pcse.decorators import prepare_rates, prepare_states
from pcse.lintul import lintul3lib
from pcse.lintul.lintul3lib import notNull, INSW, REAAND
from numpy.ma.core import exp
from numbers import Number
from pcse.lintul.stateVariables import StateVariables
from pcse.util import astro, penman
from pcse.crop.phenology import DVS_Phenology as Phenology
from pcse.crop.partitioning import DVS_Partitioning as Partitioning


class SubModel(SimulationObject):
    
    initialValues = None
    onOutput = None
    output = {}
    headerPrinted = False
    
    OUTPUT_VARS = ["TIME"] + sorted( ["WAI", "DVS", "TSUM", "TAGBM", "WST", "WLVG", "WLVD", "WSO", "LAI", "NTAC", "WRT", 
               "GTSUM", "CBALAN", "TRANRF", "NNI", "SLA", "FRACT", "FRTWET", "FLVT", "FSTT", "FSOT", 
               "RWLVG", "RWST", "RWRT", "RWSO", "CUMPAR", "LUECAL", "NUPTT", "TTRAN", "TEVAP", "PEVAP", 
               "NBALAN", "WATBAL", "NUPTR", "TNSOIL", "NDEMTO", "RNSOIL", "FERTN", "FERTNS", "WA", 
               "TIRRIG", "TRAIN", "TEXPLO", "TRUNOF", "TDRAIN"])

#     OUTPUT_VARS = ["TIME"] + sorted( ["CHECK", "FRTWET", "FLVT","FSTT","FSOT"])
     
    @classmethod
    def doOutput(cls, model, time, localVariables):
        '''
        outputs output
        :param time:
        :param localVariables: locals().copy()
        '''
        if SubModel.onOutput != None:
            
            # delegate output printing:
            if SubModel.output.has_key("TIME") and (SubModel.output["TIME"] != time ):
                rcd = []
                for v in SubModel.OUTPUT_VARS:
                    try:
                        value = SubModel.output[v]
                    except:
                        value = "-"
                    rcd.append(value)
                    
                SubModel.onOutput(rcd, SubModel.OUTPUT_VARS if not SubModel.headerPrinted else None)
                SubModel.headerPrinted = True
                SubModel.output        = {}
                 
                
            # compile value record:              
            SubModel.output["TIME"] = time
            
            for v in SubModel.OUTPUT_VARS:
                if localVariables.has_key(v):
                    SubModel.output[v] = localVariables[v]
                elif hasattr(model.states, v):
                    SubModel.output[v] = getattr(model.states, v)
                
            
        
        

class Lintul3(SubModel):
    """
* ORIGINAL COPYRGIGHT NOTICE:    
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
        DVSDR = Float(-99.)
        DVSNLT= Float(-99.)
        DVSNT = Float(-99.)
        FNTRT = Float(-99.)
        FRNX  = Float(-99.)
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
        WCAD  = Float(-99.)
        WCFC  = Float(-99.)
        WCI   = Float(-99.)
        WCST  = Float(-99.)
        WCWET = Float(-99.)
        WCWP  = Float(-99.)
        WMFAC = Float(-99.)
        RDRNS = Float(-99.)
        RDRRT = Float(-99.)
        

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
        


        
    class Lintul3States(StateVariables):
        LAI   = Float(-99.)
        ANLV  = Float(-99.)
        ANST  = Float(-99.)
        ANRT  = Float(-99.)
        ANSO  = Float(-99.)
        NUPTT = Float(-99.)
        NLOSSL= Float(-99.)
        NLOSSR= Float(-99.)
        WLVG  = Float(-99.)
        WLVD  = Float(-99.)
        WST   = Float(-99.)
        WSO   = Float(-99.)
        WRT   = Float(-99.)
        ROOTD = Float(-99.)
        GTSUM = Float(-99.) # ?? obsolete ??
        WDRT  = Float(-99.)
        CUMPAR= Float(-99.)


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
            
        
        
    # predefined attributes:        
    DSLR  = 0 
    FSHMOD = 0.0       
    
    # sub-model components for crop simulation
    pheno = Instance(SimulationObject)
    part  = Instance(SimulationObject)
    

    
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
        initialStates       = self.Lintul3States.initialValues() 

        
        # Initialize state variables
#         initialStates["TSUM"]   = init.TSUMI
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

        # Initialize components of the crop
        self.pheno = Phenology(day, kiosk, parvalues)
        self.part  = Partitioning(day, kiosk, parvalues)


        self.states = self.Lintul3States(kiosk, publish=["LAI", "ROOTD"], **initialStates)                
        self.states.initialize()
        
        kiosk.register_variable(self, "EMERG", type="R", publish=True)
        kiosk.register_variable(self, "EVAP", type="R", publish=True)
        kiosk.register_variable(self, "TRAN", type="R", publish=True)
        
        
                
    def find_subroutines(self):
        self.deathl = lintul3lib.DEATHL
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
        s = self.states
        i = self.initialValues
        
        TIME =  day.timetuple().tm_yday
        DELT = 1 # ???
        
        # Variables supplied by the weather system
        LAT              = drv.LAT
        RDD              = drv.IRRAD # ???
        TMMN             = drv.TMIN
        TMMX             = drv.TMAX
        VP               = drv.VAP  / 10 # hPa --> kPa correction from Cabo Weather
        WN               = drv.WIND
        RAIN             = drv.RAIN * 10 # cm  --> mm CORRECTION FOR NON-STANDARD cm in CABO-WEATHER

        DTR    = RDD/1.E+6  # Actual daily total global radiation (DTR, J m-2 d-1, the factor 1.E06 converts J into MJ)
        DAVTMP = 0.5 * (TMMN + TMMX)
        
        # Phenology
        self.pheno.calc_rates(day, drv)
        crop_stage = self.pheno.get_variable("STAGE")

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if crop_stage == "emerging":
            return
                

        # **********************Calling Subroutines******************************
        #  Calling the subroutine for calculating the astrological daylength.
        DAYL = astro(day, LAT, drv.IRRAD).DAYL
        
        # Calling the subroutine for converting TSUM to the developmental stage.
        DVS = self.pheno.get_variable("DVS") # self.kiosk["DVS"]
        
        # Calling the subroutine for translocatable N in leaves, stem, roots and
        # storage organs.
        ATNLV, ATNST, ATNRT, ATN = self.ntrloc(s.ANLV, s.ANST, s.ANRT, s.WLVG, s.WST, s.WRT, 
                                               p.RNFLV, p.RNFST, p.RNFRT, p.FNTRT)
        
        
                # Calling the subroutine for relative growth rate of roots, leaves, stem
        # and storage organs.
        RWLVG, RWRT, RWST, RWSO = self.relgr(TIME, p.DOYEM, EMERG, 
                                             i.WLVGI, i.WRTLI, i.WSTI, i.WSOI, 
                                             GTOTAL, FLV, FRT, FST, FSO, DLV, DRRT, DELT)
        
# Leaf growth and LAI.
        GLV    = FLV * GTOTAL
        
        # Calling the subroutine for daily increase of leaf area index.
        GLAI = self.gla(TIME, p.DOYEM, DTEFF, i.LAII, p.RGRL, DELT, SLA, s.LAI, GLV, p.NLAI, WC, p.WCWP, DVS, TRANRF, NNI)
        
        # Net rate of change of Leaf area.
        RLAI   = GLAI - DLAI
        
        # Relative death rate of leaves due to senescence/ageing.
        RDRTMP = p.RDRT(DAVTMP)
        
        # *Total vegetative biomass.
        TBGMR = s.WLVG + s.WST
        
        # Total N in green matter of the plant.
        NUPGMR = s.ANLV + s.ANST
        
        # Total leaf weight.
        WLV    = s.WLVG + s.WLVD
        
        #  Water content in the rootzone
        WA = self.kiosk["WA"]
        WC  = 0.001* WA /notNull(s.ROOTD)
        
        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        NMAXLV = p.NMXLV(DVS)
        
        # * Total above ground biomass
        TAGBM = WLV + s.WST + s.WSO ## --> output @UnusedVariable
        
        EMERG  = REAAND(TIME - p.DOYEM + 1., WC - p.WCWP)* INSW(-s.LAI, 1., 0.)
        
        # Calling the subroutine for Potential evaporation and transpiration.
        # RLWN, NRADC, PENMRC, PENMD, ...
        PEVAP, PTRAN = self.penman(DAVTMP, VP, DTR, s.LAI, WN)
        E0, ES0, ET0 = penman(day, LAT, drv.ELEV, drv.TMIN, drv.TMAX, drv.IRRAD, drv.VAP, drv.WIND, ANGSTA=-0.18, ANGSTB= -0.55)
        # NB: idem in [cm] uit drv te halen....
        
        # *Average residual N concentration.
        NRMR   = (s.WLVG * p.RNFLV + s.WST * p.RNFST) / notNull(TBGMR)
        
        #  Nitrogen uptake limiting factor at low moisture conditions in the
        #  rooted soil layer before anthesis. After anthesis there is no
        #  uptake from the soil anymore.
        NLIMIT = INSW(DVS - p.DVSNLT, INSW(WC - p.WCWP, 0., 1.), 0.0)
        
        # -------- Growth rates and dry matter production of plant organs-------*
        #  Biomass partitioning functions under (water and nitrogen)non-stressed
        #  situations
        pf = self.part.calc_rates(day, drv)
        FRTWET = pf.FR
        FLVT   = pf.FL
        FSTT   = pf.FS
        FSOT   = pf.FO
         
        # * For calculation of LUE (it can be removed once after checking, however, 
        #   I have put these new additions in these two (c and N) balances sections
        PAR    = DTR * 0.50
        DTEFF  = max ( 0., DAVTMP - p.TBASE )
#         RTSUM  = DTEFF * EMERG
        
        # Calling the subroutine for actual rates of evaporation and
        # transpiration.
        EVAP, TRAN, self.DSLR = self.evaptr(PEVAP, PTRAN, s.ROOTD, WA, 
                                            p.WCAD, p.WCWP, p.WCFC, p.WCWET, p.WCST, 
                                            p.TRANCO, DELT, p.WMFAC, RAIN, self.DSLR)
        NMAXST = p.LSNR * NMAXLV
        NMAXRT = p.LRNR * NMAXLV
        
        RROOTD = min(p.RRDMAX * INSW(WC - p.WCWP, 0., 1.) * EMERG,  p.ROOTDM - s.ROOTD)
        
        # Growth reduction function for water stress(actual trans/potential)
        TRANRF = TRAN / notNull(PTRAN)
        
        # Calling the subroutine for maximum nitrogen concentration of leaves
        # and stem.
        NOPTLV, NOPTST = self.noptm(p.FRNX, NMAXLV, NMAXST)
        
        # Maximum N content in the plant.
        NOPTS = NOPTST * s.WST
        NOPTL = NOPTLV * s.WLVG
        NOPTMR = (NOPTL+ NOPTS)/notNull(TBGMR)
        
        # Calling the subroutine for Nitrogen Nutrition Index.
        NFGMR  = NUPGMR / notNull(TBGMR)
        NNI = self.nnindx(TIME, p.DOYEM, EMERG, NFGMR, NRMR, NOPTMR)
        
        # Calling the subroutine for relative modification for root and shoot
        # allocation.
        self.FSHMOD, FRT, FLV, FST, FSO = self.subpar(p.NPART, TRANRF, NNI, FRTWET, FLVT, FSTT, FSOT, self.FSHMOD)
        
        # Calling the subroutine for total growth rate.
        # PARINT, ... 
        GTOTAL = self.growth(TIME, p.DOYEM, DTR, p.K, p.NLUE, s.LAI, p.LUE, TRANRF, NNI)
        
        # -------------Functions and parameters for rice----------------------
        
        # Specific Leaf area(m2/g).
        SLA = p.SLAC * p.SLACF(DVS) * exp(-p.NSLA * (1.-NNI))
        
        # Calling the subroutine for relative death rate of leaves.
        TSUM = self.pheno.get_variable("TSUM")
        DLV, DLAI = self.deathl(TIME, p.DOYEM, TSUM, p.TSUMAG, RDRTMP, 
                                p.RDRSHM, s.LAI, p.LAICR, s.WLVG, p.RDRNS, NNI, SLA)
        
        # Calling the subroutine for N loss due to death of leaves and roots.
        DRRT, RNLDLV, RNLDRT = self.rnld(DVS, s.WRT, p.RDRRT, p.RNFLV, DLV, p.RNFRT, p.DVSDR)
        
        # Calling the subroutine for N demand of leaves, roots and stem storage
        # organs.
        NDEML, NDEMS, NDEMR, NDEMSO = self.ndemnd(NMAXLV, NMAXST, NMAXRT, p.NMAXSO, 
                                                  s.WLVG, s.WST, s.WRT, s.WSO, RWLVG, RWST, RWRT, RWSO, 
                                                  s.ANLV, s.ANST, s.ANRT, s.ANSO, p.TCNT, DELT)
        
        #     Total Nitrogen demand.
        NDEMTO = max(0.0, (NDEML + NDEMS + NDEMR))
        
        # N supply to the storage organs.
        NSUPSO = INSW (DVS - p.DVSNT, 0., ATN / p.TCNT)

        # Rate of N uptake in grains.
        RNSO =  min(NDEMSO, NSUPSO)
        
        # --------------------Fertilizer N---------------------------------------
        #     Total nutrient uptake.
        TNSOIL = self.kiosk["TNSOIL"]
        NUPTR = (max(0., min(NDEMTO, TNSOIL))* NLIMIT ) / DELT * INSW(TIME - p.DOYEM, 0., 1.)
        
        # Calling the subroutine for N translocated from leaves, stem, and roots.
        RNTLV, RNTST, RNTRT = self.ntrans(RNSO, ATNLV, ATNST, ATNRT, ATN)
        
        # Calling the subroutine to compute the partitioning of the total
        # N uptake rate (NUPTR) over the leaves, stem and roots.
        RNULV, RNUST, RNURT = self.rnusub(TIME, p.DOYEM, EMERG, NDEML, NDEMS, NDEMR, NUPTR, NDEMTO, 
                                          i.ANLVI, i.ANSTI, i.ANRTI, DELT)
        
        RNST = RNUST-RNTST
        RNRT = RNURT-RNTRT-RNLDRT
        
        # *Rate of change of N in organs
        RNLV = RNULV-RNTLV-RNLDLV


        # *  Carbon, Nitrogen, Water balance
        CBALAN = (s.GTSUM + (i.WRTLI + i.WLVGI + i.WSTI + i.WSOI)     # @UnusedVariable
                         - (WLV + s.WST + s.WSO + s.WRT + s.WDRT))   

        NBALAN = (s.NUPTT + (i.ANLVI + i.ANSTI + i.ANRTI + i.ANSOI)   # @UnusedVariable
                         - (s.ANLV + s.ANST + s.ANRT + s.ANSO + s.NLOSSL + s.NLOSSR)) 

        self.kiosk.set_variable(self, "EMERG", EMERG)
        self.kiosk.set_variable(self, "EVAP", EVAP)
        self.kiosk.set_variable(self, "TRAN", TRAN)
        self.kiosk.set_variable(self, "NLIMIT", NLIMIT)
        self.kiosk.set_variable(self, "NUPTR", NUPTR)
        
        s.rLAI    = RLAI  
        s.rANLV   = RNLV  
        s.rANST   = RNST  
        s.rANRT   = RNRT  
        s.rANSO   = RNSO  
        s.rNUPTT  = NUPTR 
        s.rNLOSSL = RNLDLV
        s.rNLOSSR = RNLDRT
        s.rWLVG   = RWLVG 
        s.rWLVD   = DLV   
        s.rWST    = RWST  
        s.rWSO    = RWSO  
        s.rWRT    = RWRT  
        s.rROOTD  = RROOTD
        s.rGTSUM  = GTOTAL
        s.rWDRT   = DRRT  
        s.rCUMPAR = PAR   
        
        self.doOutput(self, TIME, locals().copy())


    
    
    @prepare_states
    def integrate(self, day):
        states = self.states
        delta = 1.
        
        # if before emergence there is no need to continue
        # because only the phenology is running.
        # Just run a touch() to to ensure that all state variables are available
        # in the kiosk
        # crop stage before integration
        crop_stage = self.pheno.get_variable("STAGE")
        self.pheno.integrate(day)
        if crop_stage == "emerging":
            self.touch()
            return

        # Partitioning
        self.part.integrate(day)
                        
        for s in states.listIntegratedStates():
            rate = getattr(states, 'r' + s)
            state = getattr(states, s)
            newvalue = state + delta * rate
            setattr(states, s, newvalue)

        
        
    @prepare_rates
    def updateRates(self):
        self.rates.rTSUM = 1
        
        
        
if (__name__ == "__main__"):
    from start import Lintul3Model
    
    class P:            
        __lineBuffer = {}
        __headerBuffer = {}
        __headerPrinted = False


        def printRow(self, values):
            for v in values:
                if isinstance(v, Number):
                    print "%f\t" % (v),
                else:
                    print "%s\t" % (v),
            print
            
            
        def __call__(self, values, header = None):
            if header:
                self.printRow(header)
            self.printRow(values)


        
    p = P()
    sim = Lintul3Model.start(year=1987, outputProc=p)
    l = sim.crop
    
    sim.run(365)
    SubModel.doOutput(l, 999, [])
    print "END STOP ENDJOB "*5
