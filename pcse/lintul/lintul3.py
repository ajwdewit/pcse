# -*- coding: utf-8 -*-
# Copyright (c) 2004-2015 Alterra, Wageningen-UR
# Wim de Winter(Wim.deWinter@wur.nl), April 2015
'''
LINTUL3
'''
from pcse.base_classes import SimulationObject, ParamTemplate
from pcse.base_classes import StatesWithImplicitRatesTemplate as StateVariables
from pcse.traitlets import Float, AfgenTrait, Instance
from pcse.decorators import prepare_rates, prepare_states
from pcse.util import limit
from pcse.crop.phenology import DVS_Phenology as Phenology
from pcse.lintul.lintul3lib import notNull, INSW, REAAND
from numpy.ma.core import exp, sqrt
from numbers import Number


class SubModel(SimulationObject):
    """
    Super class to enable simple output
    """
    
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

LINTUL3 is a crop model that calculates biomass production based on intercepted photosynthetically
active radiation (PAR) and light use efficiency (LUE). It is an adapted version of LINTUL2 (that simulates
potential and water-limited crop growth), including nitrogen limitation. Nitrogen stress in the model is
defined through the nitrogen nutrition index (NNI): the ratio of actual nitrogen concentration and critical
nitrogen concentration in the plant. The effect of nitrogen stress on crop growth is tested in the model
either through a reduction in LUE or leaf area (LA) or a combination of these two and further evaluated
with independent datasets. However, water limitation is not considered in the present study as the
crop is paddy rice. This paper describes the model for the case of rice, test the hypotheses of N stress
on crop growth and details of model calibration and testing using independent data sets of nitrogen
treatments (with fertilizer rates of 0 - 400 kgNha-1) under varying environmental conditions in Asia.
Results of calibration and testing are compared graphically, through Root Mean Square Deviation (RMSD),
and by Average Absolute Deviation (AAD). Overall average absolute deviation values for calibration and
testing of total aboveground biomass show less than 26% mean deviation from the observations though
the values for individual experiments show a higher deviation up to 41%. In general, the model responded
well to nitrogen stress in all the treatments without fertilizer application as observed, but between
fertilized treatments the response was varying.

Nitrogen demand, uptake and stress 

At sub-optimal nitrogen availability in the soil, nitrogen demand of the crop 
cannot be satisfied, which leads to sub-optimal crop nitrogen concentration. The 
crop nitrogen concentration below which a crop experiences nitrogen stress is 
called the critical nitrogen concentration. Nitrogen stress results in reduced 
rates of biomass production and eventually in reduced yields. Actual N content 
is the accumulated N above residual (which forms part of the cell structure). 
The critical N content is the one corresponding to half of the maximum. Nitrogen 
contents of these three reference points include those of leaves and stems, 
whereas roots are not considered since N contents of above-ground (green) parts 
are more important for photosynthesis, because of their chlorophyll content. 
However, calculation of N demand and N uptake also includes the belowground 
part. 

http://models.pps.wur.nl/sites/models.pps.wur.nl/files/LINTUL-N-Shibu-article.pdf
"""
    # Parameters, rates and states which are relevant at the main crop simulation level
    class Parameters(ParamTemplate):
        DOYEM  = Float(-99)    # The day of the year on which crop emerges.
        DVSDR  = Float(-99)    # Development stage above which deathOfLeaves of leaves and roots start.
        DVSNLT = Float(-99)    # development stage N-limit
        DVSNT  = Float(-99)    # development stage N-threshold
        FNTRT  = Float(-99)    # Nitrogen translocation from roots as a fraction of the total amount of nitrogen translocated from leaves and stem.
        FRNX   = Float(-99)    # Optimal N concentration as the fraction of maximum N concentration.
        K      = Float(-99)    # light extinction coefficient
        LAICR  = Float(-99)    # (oC d)-1, critical LAI above which mutual shading of leaves occurs,
        LRNR   = Float(-99)    # 
        LSNR   = Float(-99)    # 
        LUE    = Float(-99)    # Light use efficiency.
        NFRLVI = Float(-99)    # Initial fraction of N (g N g-1 DM) in leaves.
        NFRRTI = Float(-99)    # Initial fraction of N (g N g-1 DM) in roots.
        NFRSTI = Float(-99)    # Initial fraction of N (g N g-1 DM) in stem.
        NLAI   = Float(-99)    # Coefficient for the effect of N stress on LAI reduction(during juvenile phase)
        NLUE   = Float(-99)    # Extinction coefficient for  Nitrogen distribution down the canopy
        NMAXSO = Float(-99)    # 
        NPART  = Float(-99)    # Coefficient for the effect of N stress on leaf biomass reduction
        NSLA   = Float(-99)    # Coefficient for the effect of N stress on SLA reduction
        RDRSHM = Float(-99)    # and the maximum relative deathOfLeaves rate of leaves due to shading.
        RGRL   = Float(-99)    # Relative totalGrowthRate rate of LAI at the exponential totalGrowthRate phase
        RNFLV  = Float(-99)    # Residual N concentration in leaves
        RNFRT  = Float(-99)    # Residual N concentration in roots.
        RNFST  = Float(-99)    # Residual N concentration in stem
        ROOTDM = Float(-99)    # Maximum root depth for a rice crop.
        RRDMAX = Float(-99)    # Maximum rate of increase in rooting depth (m d-1) for a rice crop.
        SLAC   = Float(-99)    # Specific leaf area constant.
        TBASE  = Float(-99)    # Base temperature for spring wheat crop.
        TCNT   = Float(-99)    # Time coefficient(days) for N translocation.
        TRANCO = Float(-99)    # Transpiration constant (mm/day) indicating the level of drought tolerance of the wheat crop.
        TSUMAG = Float(-99)    # Temperature sum for ageing of leaves
        TSUMAN = Float(-99)    # Temperature sum for anthesis  [corresponds to TSUM1= Float(-99.)# Temperature sum emergence to anthesis]
        WCAD   = Float(-99)    # Soil hydraulic properties
        WCFC   = Float(-99)    # Soil hydraulic properties
        WCI    = Float(-99)    # Initial water content in cm3 of water/(cm3 of soil).
        WCST   = Float(-99)    # Soil hydraulic properties
        WCWET  = Float(-99)    # Soil hydraulic properties
        WCWP   = Float(-99)    # Soil hydraulic properties
        WMFAC  = Float(-99)    # water management (0                                =irrigated up to the field capacity, 1= irrigated up to saturation)
        RDRNS  = Float(-99)    # Relative deathOfLeaves rate of leaves due to N stress.
        RDRRT  = Float(-99)    # Relative deathOfLeaves rate of roots.
        RDRRT  = Float(-99)    # Relative deathOfLeaves rate of roots.        

        FLVTB  = AfgenTrait()  # Partitioning coefficients
        FRTTB  = AfgenTrait()  # Partitioning coefficients
        FSOTB  = AfgenTrait()  # Partitioning coefficients
        FSTTB  = AfgenTrait()  # Partitioning coefficients
        NMXLV  = AfgenTrait()  # Maximum N concentration in the leaves as a function of development stage.
        PHOTTB = AfgenTrait()  # Function to include the effect of photoperiodicity
        RDRT   = AfgenTrait()  # 
        SLACF  = AfgenTrait()  # Leaf area correction function as a function of development stage, DVS.        
        
        FERTAB = AfgenTrait()  # Fertilizer application as a function of TIME (g N m-2).
        NRFTAB = AfgenTrait()  # Fertilizer nitrogen recovery fraction
        ROOTDI  = Float(-99)   # initial rooting depth [m] 


        
    class Lintul3States(StateVariables):
        LAI   = Float(-99.) # leaf area index
        ANLV  = Float(-99.) # Actual N content in leaves
        ANST  = Float(-99.) # Actual N content in stem 
        ANRT  = Float(-99.) # Actual N content in root
        ANSO  = Float(-99.) # Actual N content in storage organs
        NUPTT = Float(-99.) # Total uptake of N over time (g N m-2)
        NLOSSL= Float(-99.) # total N loss by leaves
        NLOSSR= Float(-99.) # total N loss by roots
        WLVG  = Float(-99.) # Weight of green leaves
        WLVD  = Float(-99.) # Weight of dead leaves
        WST   = Float(-99.) # Weight of stem
        WSO   = Float(-99.) # Weight of storage organs 
        WRT   = Float(-99.) # Weight of roots
        ROOTD = Float(-99.)
        GTSUM = Float(-99.) # ?? obsolete ??
        WDRT  = Float(-99.) # dead roots (?)
        CUMPAR= Float(-99.)
        TNSOIL= Float(-99.) # Amount of inorganic N available for crop uptake.



    class InitialValues(object):
        
        def __init__(self, parameters):
            # Read initial states
            self.ROOTDI= parameters.ROOTDI
            self.TSUMI = 0.0
            self.WLVGI = 2.4
            self.WRTLI = 3.6
            self.WSOI  = 0.0
            self.WSTI  = 0.0 
            self.WCI   = 0.40   # Initial water content in cm3 of water/(cm3 of soil).   
            
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
    
    # sub-model components for crop simulation
    pheno = Instance(SimulationObject)
    #     part  = Instance(SimulationObject)
    

    
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
        #         self.part  = Partitioning(day, kiosk, parvalues)


        self.states = self.Lintul3States(kiosk, publish=["LAI", "ROOTD"], **initialStates)                
        self.states.initialize()
        
        kiosk.register_variable(self, "EVAP",  type="R", publish=True)
        kiosk.register_variable(self, "TRAN",  type="R", publish=True)
        kiosk.register_variable(self, "RROOTD",  type="R", publish=True)
        
        
                
    @prepare_rates
    def calc_rates(self, day, drv):
        # dynamic calculations
        p = self.params
        s = self.states
        i = self.initialValues
        
        DELT    = 1 # ???

        WA      = self.kiosk["WA"]
        DVS     = self.pheno.get_variable("DVS")    # self.kiosk["DVS"]
        TSUM    = self.pheno.get_variable("TSUM")
        
        # Variables supplied by the weather system
        RDD     = drv.IRRAD # ???
        TMMN    = drv.TMIN
        TMMX    = drv.TMAX
        RAIN    = drv.RAIN * 10 # cm  --> mm CORRECTION FOR NON-STANDARD cm in CABO-WEATHER

        DTR     = RDD/1.E+6  # Actual daily total global radiation (DTR, J m-2 d-1, the factor 1.E06 converts J into MJ)
        PAR     = DTR * 0.50
        DAVTMP  = 0.5 * (TMMN + TMMX)
        DTEFF   = max ( 0., DAVTMP - p.TBASE )
        
        #  calculating the astrological daylength.
        """
        A new sub-routine has been added in LINTUL3 to calculate the daylength, which 
        affects crop development. Normally, photoperiodic daylength exceeds astronomical 
        daylength (Wormer, 1954; Vergara and Chang, 1985). Photoperiodic daylength is 
        calculated in the model as a function of solar elevation (angle of sun above horizon), 
        determined by latitude and day of the year (Goudriaan and Van Laar, 1994). 
        Photoperiodsensitivity of the rice crop, defined as a function of daylength is 
        included in the model by modifying the daily increment of the heat sum. 
        """
        #         DAYL    = astro(day, drv.LAT, drv.IRRAD).DAYL # obsolete
        
        # potential rates of evaporation and transpiration:
        PEVAP, PTRAN = self.potentialEvapoTranspiration(drv)
               
        # actual rates of evaporation and transpiration:
        EVAP, TRAN   = self.evapotranspiration(PEVAP, PTRAN, WA, RAIN, DELT) 

        # ----------Emergence,Temperature sum and Developmental stages----------*
        # Water content in the rootzone
        WC      = 0.001* WA /notNull(s.ROOTD)

        """
        Crop phenology
        
        Crop development, i.e. the order and rate of appearance of vegetative and 
        reproductive organs, is defined in terms of phenological developmental stage 
        (DVS) as a function of heat sum, which is the cumulative daily effective 
        temperature. Daily effective temperature is the average temperature above a 
        crop-specific base temperature (for rice 8C). Some rice varieties are 
        photoperiodsensitive, i.e. flowering depends on the length of the light period 
        during the day in addition to the temperature during the vegetative stage. 
        """
        self.pheno.calc_rates(day, drv)
        crop_stage = self.pheno.get_variable("STAGE")

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if crop_stage == "emerging":
            return
                
        # code below is executed only POST-emergence
            
        # translocatable N in leaves, stem, roots and
        # storage organs.
        ATNLV, ATNST, ATNRT, ATN = self.translocatable_N()
        
        # Relative deathOfLeaves rate of leaves due to senescence/ageing.
        RDRTMP  = p.RDRT(DAVTMP)
        
        # *Total vegetative biomass.
        TBGMR   = s.WLVG + s.WST
        
        # *Average residual N concentration.
        NRMR    = (s.WLVG * p.RNFLV + s.WST * p.RNFST) / notNull(TBGMR)
        
        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        NMAXLV  = p.NMXLV(DVS)        
        NMAXST  = p.LSNR * NMAXLV
        NMAXRT  = p.LRNR * NMAXLV
        
        # maximum nitrogen concentration of leaves and stem.
        NOPTLV  = p.FRNX * NMAXLV
        NOPTST  = p.FRNX * NMAXST
                
        # Maximum N content in the plant.
        NOPTS   = NOPTST * s.WST
        NOPTL   = NOPTLV * s.WLVG
        NOPTMR  = (NOPTL+ NOPTS)/notNull(TBGMR)
        
        # Total N in green matter of the plant.
        NUPGMR  = s.ANLV + s.ANST
        
        # Nitrogen Nutrition Index.
        """
        Nitrogen stress 

        A crop is assumed to experience N stress at N concentrations below a critical 
        value for unrestricted growth. To quantify crop response to nitrogen shortage, a 
        Nitrogen Nutrition Index (NNI) is defined, ranging from 0 (maximum N shortage) 
        to 1 (no N shortage): 
        
        NNI = (actual crop [N] - residual [N]) / (critical [N] - residual [N]) 
        
        Critical crop nitrogen concentration, the lower limit of canopy nitrogen 
        concentration in leaves and stems required for unrestricted growth, has been 
        taken as half the maximum nitrogen concentration. An experimental basis for such 
        an assumption can be derived from the study of Zhen and Leigh (1990), who 
        reported that nitrate accumulation in plant occurs in significant quantity when 
        the N needs to reach the maximum growth were fulfilled and the mean value of 
        nitrate accumulated beyond the criticalNconcentration was about 50% for 
        different stages. 
        """
        NFGMR   = NUPGMR / notNull(TBGMR)
        NNI     = limit (0.001, 1.0, ((NFGMR-NRMR)/notNull(NOPTMR-NRMR))) 

        # -------- Growth rates and dry matter production of plant organs-------*
        #  Biomass partitioning functions under (water and nitrogen)non-stressed
        #  situations
        """
        Biomass partitioning
        
        Biomass formed at any time during crop growth is partitioned amongits organs 
        (Fig. 1), i.e. roots, stems, leaves and storage organs, with partitioning 
        factors defined as a function of development stage (Fig. 2) (Drenth et al., 
        1994), which thus provides the rates of growth of these organs: 
        
        dW/dt[i] = Pc[i] * dW / dt 
        
        where (dW/dt) is the rate of biomass growth (gm-2 d-1); (dW/dt)[i] and Pc[i] are 
        the rate of growth (gm-2 d-1) of and the biomass partitioning factor to organ i 
        (g organ-i g-1 biomass), respectively. Leaf, stem and root weights of the 
        seedlings at the time of transplanting are input parameters for the model. The 
        time course of weights of these organs follows from integration of their net 
        growth rates, i.e. growth rates minus death rates, the latter being defined as a 
        function of physiological age, shading and stress. 
        """
        FRTWET  = p.FRTTB( DVS )
        FLVT    = p.FLVTB( DVS )
        FSTT    = p.FSTTB( DVS)
        FSOT    = p.FSOTB( DVS )        
        
        
        """
        Leaf area development
        
        The time course of LAI is divided into two stages: an exponential stage during 
        the juvenile phase, where leaf area development is a function of temperature, 
        and a linear stage where it depends on the increase in leaf biomass (Spitters, 
        1990; Spitters and Schapendonk, 1990). The death of leaves due to senescence 
        that may be enhanced by shading and/or stress leads to a corresponding loss in 
        leaf area. The specific leaf area is used for the conversion of dead leaf 
        biomass to corresponding loss in leaf area. The death of leaves due to 
        senescence occurs only after flowering and the rate depends on crop age 
        (function adopted from ORYZA2000, Bouman et al., 2001). The excessive growth of 
        leaves also result in death of leaves due to mutual shading. The death of leaves 
        due to shading is determined by a maximum death rate and the relative proportion 
        of leaf area above the critical LAI (4.0) (Spitters, 1990; Spitters and 
        Schapendonk, 1990). The net rate of change in leaf area (dLAI/dt) is the 
        difference between its growth rate and its death rate: 
        
        dLAI/dt = dGLAI / dt - dDLAI / dt

        where (dGLAI/dt) is the leaf area growth rate and (dDLAI/dt) is the
        leaf area death rate.
        """
        # Specific Leaf area(m2/g).
        SLA = p.SLAC * p.SLACF(DVS) * exp(-p.NSLA * (1.-NNI))
        
        # Growth reduction function for water stress(actual trans/potential)
        TRANRF = TRAN / notNull(PTRAN)
        
        # relative modification for root and shoot allocation.
        FRT, FLV, FST, FSO = self.dryMatterPartitioningFractions(p.NPART, TRANRF, NNI, FRTWET, FLVT, FSTT, FSOT)
        
        # total totalGrowthRate rate.
        GTOTAL = self.totalGrowthRate(DTR, TRANRF, NNI)

        # Leaf totalGrowthRate and LAI.
        GLV    = FLV * GTOTAL
        
        # daily increase of leaf area index.
        GLAI = self.gla(DTEFF, i.LAII, DELT, SLA, GLV, WC, DVS, TRANRF, NNI)
        
        # relative deathOfLeaves rate of leaves.
        DLV, DLAI = self.deathRateOfLeaves(TSUM, RDRTMP, NNI, SLA)
        
        # Net rate of change of Leaf area.
        RLAI   = GLAI - DLAI
        
        # Root totalGrowthRate
        """
        Root growth
        
        The root system is characterized by its vertical extension in the soil profile. 
        At emergence or at transplanting for transplanted rice. Effect of N stress (NNI) 
        on crop growth through its effect on LA and LUE. rooting depth is initialized. 
        Roots elongate at a constant daily rate, until flowering, provided soil water 
        content is above permanent wilting point (PWP), whereas growth ceases when soil 
        is drier than PWP or when a certain preset maximum rooting depth is reached 
        (Spitters and Schapendonk, 1990; Farr� et al., 2000). However, in rice in a 
        flooded situation, soil will always be at saturation and, therefore, the maximum 
        rooting depth corresponds to the physiological maximum, which is taken as 0.7m 
        """
        RROOTD = min(p.RRDMAX,  p.ROOTDM - s.ROOTD) if (WC > p.WCWP) else 0.0
        
        # N loss due to deathOfLeaves of leaves and roots.
        DRRT    = 0. if (DVS < p.DVSDR) else s.WRT * p.RDRRT        
        RNLDLV  = p.RNFLV* DLV
        RNLDRT  = p.RNFRT* DRRT
        
        # relative totalGrowthRate rate of roots, leaves, stem
        # and storage organs.
        RWLVG, RWRT, RWST, RWSO = self.relativeGrowthRates(GTOTAL, FLV, FRT, FST, FSO, DLV, DRRT)
        
        
        """
        Nitrogen demand 

        Total crop nitrogen demand equals the sum of the nitrogen demands of its 
        individual organs (excluding storage organs, for which nitrogen demand is met by 
        translocation from the other organs, i.e. roots, stems and leaves) (Fig. 3). 
        Nitrogen demand of the individual organs is calculated as the difference 
        betweenmaximum and actual organ nitrogen contents. The maximum nitrogen content 
        is defined as a function of canopy development stage (Drenth et al., 1994). 
        Total N demand (TNdem: gm-2 d-1) of the crop is: 
        
        TNdem = sum(Nmax,i - ANi / dt) 
        
        where Nmax,i is the maximum nitrogen concentration of organ i (gN/g biomass, 
        with i referring to leaves, stems and roots), Wi is the weight of organ i 
        (g biomass/m2), andANi is the actual nitrogen content of organ i (gN/m2). 
        """

        # N demand of leaves, roots and stem storage organs.
        NDEML   =  max(NMAXLV   * s.WLVG - s.ANLV, 0.)
        NDEMS   =  max(NMAXST   * s.WST  - s.ANST, 0.)
        NDEMR   =  max(NMAXRT   * s.WRT  - s.ANRT, 0.)
        NDEMSO  =  max(p.NMAXSO * s.WSO  - s.ANSO, 0.) / p.TCNT
        
        # N supply to the storage organs.
        NSUPSO  = INSW (DVS - p.DVSNT, 0., ATN / p.TCNT)

        # Rate of N uptake in grains.
        RNSO    =  min(NDEMSO, NSUPSO)
        
        # Total Nitrogen demand.
        NDEMTO  = max(0.0, (NDEML + NDEMS + NDEMR))
        
        """
        About 75�90% of the total N uptake at harvest takes place before
        anthesis and, in conditions of high soil fertility, post-anthesis N uptake
        may contribute up to 25% but would exclusively end up in the grain
        as protein. Therefore, this nitrogen would not play any role in the
        calculation of nitrogen stress that influences the biomass formation.
        Therefore, nitrogen uptake is assumed to cease at anthesis,
        as nitrogen content in the vegetative parts hardly increases after
        anthesis
        """
        
        #  Nitrogen uptake limiting factor at low moisture conditions in the
        #  rooted soil layer before anthesis. After anthesis there is no
        #  uptake from the soil anymore.
        NLIMIT  = 1.0 if (DVS < p.DVSNLT) and (WC >= p.WCWP) else 0.0
            
        NUPTR   = (max(0., min(NDEMTO, s.TNSOIL))* NLIMIT ) / DELT 
        
        # N translocated from leaves, stem, and roots.
        RNTLV   = RNSO* ATNLV/ notNull(ATN)
        RNTST   = RNSO* ATNST/ notNull(ATN)
        RNTRT   = RNSO* ATNRT/ notNull(ATN)
        
        # compute the partitioning of the total N uptake rate (NUPTR) over the leaves, stem and roots.
        RNULV, RNUST, RNURT = self.N_uptakeRates(NDEML, NDEMS, NDEMR, NUPTR, NDEMTO)
                
        RNST    = RNUST-RNTST
        RNRT    = RNURT-RNTRT-RNLDRT
        
        # Rate of change of N in organs
        RNLV    = RNULV-RNTLV-RNLDLV

        # ****************SOIL NITROGEN SUPPLY***********************************
        """
        Soil–crop nitrogen balance 

        The mineral nitrogen balance of the soil is the difference between nitrogen 
        added through mineralization and/or fertilizer, and nitrogen removed by crop 
        uptake and losses. The net rate of change of N in soil (dN/dt in gm-2 d-1) is: 
        
        dN/dt[soil]= Nmin + (FERTN * NRF) - dNU/dt 
        
        where Nmin is the nitrogen supply through mineralization and biological N 
        fixation, FERTN is the fertilizer nitrogen application rate, NRF is the 
        fertilizer nitrogen recovery fraction and dNU/dt is the rate of nitrogen uptake 
        by the crop, which is calculated as the minimum of the N supply from the soil 
        and the N demand from the crop. 
        """
        
        #  Soil N supply (g N m-2 d-1) through mineralization.
        RTMIN   = 0.10 * NLIMIT

        # ---------------Fertilizer application---------------------------------*
        TIME    = day.timetuple().tm_yday
        FERTN   = p.FERTAB(TIME)
        NRF     = p.NRFTAB(TIME)
        
        #  Change in inorganic N in soil as function of fertilizer
        #  input, soil N mineralization and crop uptake.
        FERTNS = FERTN * NRF
        RNSOIL = FERTNS/DELT -NUPTR + RTMIN
        
        # Total leaf weight.
        WLV     = s.WLVG + s.WLVD
        
        # Total above ground biomass
        TAGBM   = WLV + s.WST + s.WSO ## --> output @UnusedVariable
        
        # Carbon, Nitrogen, Water balance
        CBALAN = (s.GTSUM + (i.WRTLI + i.WLVGI + i.WSTI + i.WSOI)     # @UnusedVariable
                         - (WLV + s.WST + s.WSO + s.WRT + s.WDRT))   

        NBALAN = (s.NUPTT + (i.ANLVI + i.ANSTI + i.ANRTI + i.ANSOI)   # @UnusedVariable
                         - (s.ANLV + s.ANST + s.ANRT + s.ANSO + s.NLOSSL + s.NLOSSR)) 

        # export local variables for latyer use in other models
        self.kiosk.set_variable(self, "EVAP", EVAP)
        self.kiosk.set_variable(self, "TRAN", TRAN)
        self.kiosk.set_variable(self, "RROOTD", RROOTD)
        
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
        s.rTNSOIL = RNSOIL
        
        self.doOutput(self, day.timetuple().tm_yday, locals().copy())


    
    
    @prepare_states
    def integrate(self, day):
        delta = 1.
        
        # if before emergence there is no need to continue
        # because only the phenology is running.
        # Just run a touch() to to ensure that all state variables are available
        # in the kiosk

        self.pheno.integrate(day)
        if self.pheno.get_variable("STAGE") == "emerging":
            self.touch()
            return

        # Partitioning
        #         self.part.integrate(day)

        self.states.integrate(delta)              
        

        
    def potentialEvapoTranspiration(self, drv):
        # Potential evaporation and transpiration.
        ES0, ET0 = [10 * e for e in [drv.ES0, drv.ET0]]
        pevap = exp(-0.5 * self.states.LAI) * ES0
        pevap = max(0., pevap)
        ptran = (1. - exp(-0.5 * self.states.LAI)) * ET0
        ptran = max(0., ptran)
        return pevap, ptran

        
        
    # ---------------------------------------------------------------------*
    #  def EVAPTR                                                   *:
    #  Purpose: To compute actual rates of evaporation and transpiration.  *
    # ---------------------------------------------------------------------*                                                  *
    def evapotranspiration(self, PEVAP, PTRAN, WA, RAIN, DELT): 

        p = self.params
        
        # see also classic_waterbalance.py    
        WC   = 0.001 * WA   / notNull(self.states.ROOTD)
        WAAD = 1000. * p.WCAD * self.states.ROOTD
         
          
        if (RAIN >=0.5):
            EVS  = PEVAP
            self.DSLR = 1.
        else:
            self.DSLR = self.DSLR + 1.
            EVSMXT    = PEVAP*(sqrt (self.DSLR) - sqrt(self.DSLR - 1.))
            EVS       = min(PEVAP, EVSMXT + RAIN)
            
        WCCR = p.WCWP + max( 0.01, PTRAN/(PTRAN + p.TRANCO) * (p.WCFC - p.WCWP))
          
        #  If the soil is flooded for rice, : the soil is assumed to be
        #  permanently saturated and there is no effect of the high water
        #  content on crop totalGrowthRate, because rice has aerenchyma.
        #  Thus FR is formulated as below:
          
        if (p.WMFAC  >= 1.0 ):
            if (WC > WCCR):
                FR = 1.
            else:
                FR = limit( 0., 1., (WC - p.WCWP) / (WCCR - p.WCWP))
        
        #  If soil is irrigated but not flooded, : soil water content is
        #  assumed to be at field capacity and the critical water content
        #  that affects crop totalGrowthRate (FR) is formulated as below:      
        else:
            if (WC > WCCR):
                FR = limit( 0., 1., (p.WCST - WC)/(p.WCST - p.WCWET))
            else:
                FR = limit( 0., 1., (WC - p.WCWP)/(WCCR - p.WCWP))
        
        TRAN = PTRAN * FR
        
        AVAILF = min( 1., ((WA-WAAD)/DELT)/notNull(EVS+TRAN) )
        
        EVAP = EVS * AVAILF
        TRAN = TRAN * AVAILF
        
        return EVAP, TRAN
    
    
    # ---------------------------------------------------------------------*
    #  def GLA                                                      *:
    #  Purpose: This subroutine computes daily increase of leaf area index *
    #        (ha leaf/ ha ground/ d).                                   *
    # ---------------------------------------------------------------------*
    def gla(self, DTEFF, LAII,  DELT, SLA, GLV, WC, DVS, TRANRF, NNI):

        p = self.params
        LAI = self.states.LAI
        
        #---- Growth during maturation stage:
        GLAI = SLA * GLV
        
        #---- Growth during juvenile stage:
        if ((DVS  <  0.2) and (LAI  <  0.75)):
            GLAI = (LAI * (exp(p.RGRL * DTEFF * DELT) - 1.)/ DELT )* TRANRF* exp(-p.NLAI* (1.0 - NNI))
        
        #---- Growth at day of seedling emergence:
        if ((LAI == 0.) and (WC > p.WCWP)):
            GLAI = LAII / DELT  
        
        return GLAI
    
    
    #----------------------------------------------------------------------*
    #  def SUBPAR                                                  *:
    #  Purpose: Dry matter partitioning fractions: leaves, stem and storage organs.
    #----------------------------------------------------------------------*
    def dryMatterPartitioningFractions(self, NPART, TRANRF, NNI, FRTWET, FLVT, FSTT, FSOT):
      
        if(TRANRF  <  NNI):
            #  Water stress is more severe as compared to nitrogen stress and
            #  partitioning will follow the original assumptions of LINTUL2*
      
            FRTMOD = max( 1., 1./(TRANRF+0.5))
            FRT    = FRTWET * FRTMOD
            FSHMOD = (1.-FRT) / (1.-FRT/FRTMOD)
            FLV    = FLVT * FSHMOD
            FST    = FSTT * FSHMOD
            FSO    = FSOT * FSHMOD
            
        else:
            
            # Nitrogen stress is more severe as compared to water stress and the
            # less partitioning to leaves will go to the roots*
            
            FLVMOD = exp(-NPART* (1.0-NNI))
            FLV    = FLVT * FLVMOD
            MODIF  = (1.-FLV)/(1.-(FLV/FLVMOD))
            FST    = FSTT *  MODIF
            FRT    = FRTWET* MODIF
            FSO    = FSOT *  MODIF
    
        return FRT, FLV, FST, FSO # FLVMOD removed from signature - WdW
    
    
    
    # ---------------------------------------------------------------------*
    #  def GROWTH                                                   *:
    #  Purpose: To compute the total totalGrowthRate rate.                          *
    # ---------------------------------------------------------------------*
    def totalGrowthRate(self, DTR, TRANRF, NNI):
        """
        (Monteith, 1977; Gallagher and Biscoe, 1978; Monteith, 1990) have
        shown that biomass formed per unit intercepted light, LUE (Light
        Use Efficiency, g dry matter MJ-1), is relatively more stable. Hence,
        maximum daily growth rate can be defined as the product of
        intercepted PAR (photosynthetically active radiation, MJm-2 d-1)
        and LUE. Intercepted PAR depends on incident solar radiation,
        the fraction that is photosynthetically active (0.5) (Monteith and
        Unsworth, 1990; Spitters, 1990), and LAI (m2 leafm-2 soil) according
        to Lambert�Beer�s law:
        Q = 0.5Q0[1 - e(-k LAI)] 
        where Q is intercepted PAR(MJm-2 d-1), Q0 is daily global radiation
        (MJm-2 d-1), and k is the attenuation coefficient for PAR in the
        canopy.
        """
        p = self.params
        PARINT = 0.5 * DTR * (1.- exp(-p.K * self.states.LAI))
        
        GTOTAL = p.LUE * PARINT
        
        if(TRANRF  <=  NNI):
            #  Water stress is more severe as compared to nitrogen stress and
            #  partitioning will follow the original assumptions of LINTUL2*
            
            GTOTAL *= TRANRF
        
        else:
        
            #  Nitrogen stress is more severe as compared to water stress and the
            #  less partitioning to leaves will go to the roots*
            
            GTOTAL *= exp(-p.NLUE * (1.0 - NNI))
      
        return GTOTAL
    
    
    
    # ---------------------------------------------------------------------*
    #  def RELGR                                                    *:
    #  Purpose: To compute the relative totalGrowthRate rate of roots, leaves, stem *
    #        and storage organs.                                        *
    # ---------------------------------------------------------------------*
    def relativeGrowthRates(self, GTOTAL, FLV, FRT, FST, FSO, DLV, DRRT):
      
        RWLVG = GTOTAL * FLV - DLV
        RWRT  = GTOTAL * FRT - DRRT
        RWST  = GTOTAL * FST
        RWSO  = GTOTAL * FSO

        return RWLVG, RWRT, RWST, RWSO
    
    
    
    # ---------------------------------------------------------------------*
    #  def RNUSUB                                                   *:
    #  Purpose: To compute the partitioning of the total N uptake rate     *
    #        (NUPTR) over the leaves, stem, and roots.                  *
    # ---------------------------------------------------------------------*
    def N_uptakeRates(self, NDEML, NDEMS, NDEMR, NUPTR, NDEMTO):
      
        RNULV = (NDEML / notNull(NDEMTO))* NUPTR
        RNUST = (NDEMS / notNull(NDEMTO))* NUPTR
        RNURT = (NDEMR / notNull(NDEMTO))* NUPTR

        return RNULV, RNUST, RNURT
    
    
    
    # ---------------------------------------------------------------------*
    #  def NTRLOC                                                   *:
    #  Purpose: To compute the translocatable N in the organs.             *
    # ---------------------------------------------------------------------*
    def translocatable_N(self):      
        s = self.states
        p = self.params
        ATNLV = max (0., s.ANLV - s.WLVG * p.RNFLV)
        ATNST = max (0., s.ANST - s.WST  * p.RNFST)
        ATNRT = min((ATNLV + ATNST) * p.FNTRT, s.ANRT - s.WRT * p.RNFRT)
        ATN   = ATNLV +  ATNST + ATNRT
        
        return ATNLV, ATNST, ATNRT, ATN

    
    
    # ---------------------------------------------------------------------*
    #  def DEATHL                                                   *:
    #  Purpose: To compute the relative deathOfLeaves rate of leaves due to age, *
    #        shading amd due to nitrogen stress.                        *
    # ---------------------------------------------------------------------*
    def deathRateOfLeaves(self, TSUM, RDRTMP, NNI, SLA):
      
        p = self.params
        s = self.states
        
        RDRDV = 0. if (TSUM < p.TSUMAG) else RDRTMP
        
        RDRSH = max(0., p.RDRSHM * (s.LAI - p.LAICR) / p.LAICR)
        RDR   = max(RDRDV, RDRSH)
        
        if (NNI  <  1.):
            DLVNS   = s.WLVG * p.RDRNS * (1. - NNI)
            DLAINS  = DLVNS * SLA
        else:
            DLVNS   = 0.
            DLAINS  = 0.
        
        DLVS  = s.WLVG * RDR
        DLAIS = s.LAI  * RDR
        
        DLV   = DLVS + DLVNS
        DLAI  = DLAIS + DLAINS
    
        return DLV, DLAI
    






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

#     Help(Lintul3)    
    sim.run(365)
    SubModel.doOutput(l, 999, [])
    print "END STOP ENDJOB "*5
