from numpy.ma.core import sin, cos, sqrt, arcsin, exp
from pcse.util import limit



def notNull(X):
    """
    REAL FUNCTION notNull (X). This function can be used to avoid "divide
    by zero" errors in division. The function result is defined as notNull =
    X if X is not 0 and notNull = 1 when X is 0.
    """
    return X if (X <> 0.) else 1.



def INSW(X1, X2, X3):
    """
    input switch relay
    """
    return X2 if X1 < 0 else X3

       
       
def REAAND(X1, X2):
    """
    Returns 1.0 if both input values are positive, otherwise Y=0.0.
    """
    return 1. if (X1 > 0 and X2 > 0) else 0.   


    
#************************Subroutines************************************
# ---------------------------------------------------------------------*
#  def ASTRO                                                    *:
#  Purpose: This subroutine calculates the astrological daylength,     *
#        based on Goudriaan and van Laar 1994, around p.30.         *
# ---------------------------------------------------------------------*
        
def ASTRO(DOY, LAT):
    
    #-----PI VALUE
    PI     = 3.1415926
    
    # SINE AND COSINE OF LATITUDE
    SINLAT = sin (PI * LAT/180.)
    COSLAT = cos (PI * LAT/180.)
    
    # MAXIMAL SINE OF DECLINATION
    SINDCM = sin (PI * 23.45/180.)
    
    # SINE AND COSINE OF DECLINATION  (EQUATIONS 3.4, 3.5)
    
    SINDEC = -SINDCM * cos (2.* PI * (DOY+10.)/365.)
    COSDEC = sqrt (1.-SINDEC * SINDEC)
    
    # THE TERMS A AND B ACCORDING TO EQUATION 3.3
    
    A     = SINLAT * SINDEC
    B     = COSLAT * COSDEC
    
    # DAYLENGTH ACCORDING TO EQUATION 3.6
    
    DAYL  = 12. * (1. + (2./PI)* arcsin (A/B))
    return DAYL
        
# ---------------------------------------------------------------------*
#  def SUBDVS                                                   *:
#  Purpose: To compute the developmental stage corresponding to the    *
#  temperature sum.                                                    *
# ---------------------------------------------------------------------*
DVS1 = 0.
DVS2 = 0.
        
def SUBDVS(TIME, DOYEM, TSUM, TSUMAN, TSUMMT):

    global DVS1, DVS2
        
    if( TIME  <=  DOYEM ):
        DVS1 = TSUM/TSUMAN
    
    elif( TIME  >= DOYEM  and  TSUM  <=  TSUMAN ): 
        DVS1   = TSUM/TSUMAN
    
    else:
        DVS2   = (TSUM-TSUMAN)/TSUMMT
    
    DVS = DVS1 + DVS2
    return DVS


      
# ---------------------------------------------------------------------*
#  def PENMAN                                                   *:
#  Purpose: Computation of the PENMAN EQUATION.                        *
# ---------------------------------------------------------------------*
      
# DEFINE_CALL PENMAN(INPUT,INPUT,INPUT,INPUT,INPUT,   ...
#                    OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)      

def PENMAN(DAVTMP,  ## degree C 
           VP,      ## [kPa] 
           DTR,     ## [MJ] 
           LAI,     ## [-] 
           WN):     ## ?[m/s]?     
    
    DTRJM2 = DTR * 1.E6     ## [J/m2]
    BOLTZM = 5.668E-8       ## 5.66e-8 [W/m2.K4] = Stephan-Boltzmann constant
    LHVAP  = 2.4E6          ## equal to 2.4 [MJ/kg] at 30 C with only a small temperature dependence
    PSYCH  = 0.067          ## [kPa.K]
    
    BBRAD  = BOLTZM * (DAVTMP+273.)**4 * 86400.
    SVP    = 0.611 * exp(17.4 * DAVTMP / (DAVTMP + 239.))
    SLOPE  = 4158.6 * SVP / (DAVTMP + 239.)**2
    RLWN   = BBRAD * max(0., 0.55*(1.-VP/SVP))
    NRADS  = DTRJM2 * (1.-0.15) - RLWN
    NRADC  = DTRJM2 * (1.-0.25) - RLWN
    PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)
    PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)
    
    WDF    = 2.63 * (1.0 + 0.54 * WN)
    PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)
    
#     E0, ES0, ET0
#         E0      -  Penman potential evaporation from a free water surface [mm/d]
#         ES0     -  Penman potential evaporation from a moist bare soil surface [mm/d] != PENMRS
#         ET0     -  Penman potential transpiration from a crop canopy [mm/d] != PENMRC
    PEVAP  = exp(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP # (PENMRS + PENMD) / LHVAP == ES0
    PEVAP  = max (0., PEVAP)
    PTRAN  = (1.-exp(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP # (PENMRC + PENMD) / LHVAP == ET0
    PTRAN  = max( 0., PTRAN)
    
#     return RLWN, NRADC, PENMRC, PENMD, PEVAP, PTRAN
    return PEVAP, PTRAN
      
      
# ---------------------------------------------------------------------*
#  def EVAPTR                                                   *:
#  Purpose: To compute actual rates of evaporation and transpiration.  *
# ---------------------------------------------------------------------*                                                  *
      
# DEFINE_CALL EVAPTR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
#                    INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)      


def EVAPTR(PEVAP, PTRAN, ROOTD, WA, WCAD, WCWP, WCFC, WCWET, WCST,  
           TRANCO, DELT, WMFAC, RAIN, DSLR):
    # see also classic_waterbalance.py    
    WC   = 0.001 * WA   / notNull(ROOTD)
    WAAD = 1000. * WCAD * ROOTD
#     WAFC = 1000. * WCFC * ROOTD
   
    
    if (RAIN >=0.5):
        EVS  = PEVAP
        DSLR = 1.
    else:
        DSLR   = DSLR+1.
        EVSMXT = PEVAP*(sqrt (DSLR)-sqrt(DSLR-1.))
        EVS    = min (PEVAP, EVSMXT+RAIN)
      
    WCCR = WCWP + max( 0.01, PTRAN/(PTRAN+TRANCO) * (WCFC-WCWP) )
      
    #  If the soil is flooded for rice, : the soil is assumed to be
    #  permanently saturated and there is no effect of the high water
    #  content on crop growth, because rice has aerenchyma.
    #  Thus FR is formulated as below:
      
    if (WMFAC  >= 1.0 ):
        if (WC > WCCR):
            FR = 1.
        else:
            FR = limit( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )
    
    #  If soil is irrigated but not flooded, : soil water content is
    #  assumed to be at field capacity and the critical water content
    #  that affects crop growth (FR) is formulated as below:
    
    else:
        if (WC > WCCR):
            FR = limit( 0., 1., (WCST-WC)/(WCST-WCWET) )
        else:
            FR = limit( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )

    TRAN = PTRAN * FR

    AVAILF = min( 1., ((WA-WAAD)/DELT)/notNull(EVS+TRAN) )
    
    EVAP = EVS * AVAILF
    TRAN = TRAN * AVAILF

#     return WCCR, EVAP, TRAN, FR, AVAILF, DSLR
    return EVAP, TRAN, DSLR


# ---------------------------------------------------------------------*
#  def DRUNIR                                                   *:
#  Purpose: To compute rates of drainage, runoff and irrigation.       *
# ---------------------------------------------------------------------*

# DEFINE_CALL DRUNIR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
#                    INPUT,INPUT,INPUT,INPUT,  OUTPUT,OUTPUT,OUTPUT)

def DRUNIR(RAIN, EVAP, TRAN, IRRIGF,                         
            DRATE, DELT, WA, ROOTD, WCFC, WCST, WMFAC):
    
    WC   = 0.001 * WA   / notNull(ROOTD)
    WAFC = 1000. * WCFC * ROOTD
    WAST = 1000. * WCST * ROOTD
    
    DRAIN  = limit( 0., DRATE, (WA-WAFC)/DELT + (RAIN - EVAP - TRAN)                  )
    RUNOFF =          max( 0., (WA-WAST)/DELT + (RAIN - EVAP - TRAN - DRAIN)          )
    
    if (WMFAC  >= 1.0):
        #     If a soil is irrigated by flooding, : soil water content is
        #     kept at saturation via "irrigation events".
        IRRIG  = IRRIGF * max( 0., (WAST-WA)/DELT - (RAIN - EVAP - TRAN - DRAIN - RUNOFF) )
    else:
        
        #     If soil is irrigated but not flooded, : soil water content
        #     is kept at field capacity via "irrigation events".
        IRRIG  = IRRIGF * max( 0., (WAFC-WA)/DELT - (RAIN - EVAP - TRAN - DRAIN - RUNOFF) )

    return DRAIN, RUNOFF, IRRIG




# ---------------------------------------------------------------------*
#  def GLA                                                      *:
#  Purpose: This subroutine computes daily increase of leaf area index *
#        (ha leaf/ ha ground/ d).                                   *
# ---------------------------------------------------------------------*

def GLA(TIME, DOYEM, DTEFF, LAII, RGRL, DELT, SLA, LAI, GLV, NLAI,  
        WC, WCWP, DVS, TRANRF, NNI):
    
    #---- Growth during maturation stage:
    GLAI = SLA * GLV
    
    #---- Growth during juvenile stage:
    if ((DVS  <  0.2) and (LAI  <  0.75)):
        GLAI = (LAI * (exp(RGRL * DTEFF * DELT) - 1.)/ DELT )* TRANRF* exp(-NLAI* (1.0 - NNI))
    
    #---- Growth at day of seedling emergence:
    if ((TIME >= DOYEM) and (LAI == 0.) and (WC > WCWP)):
        GLAI = LAII / DELT  
    
    #---- Growth before seedling emergence:
    if (TIME < DOYEM):
        GLAI = 0.
    
    return GLAI


#----------------------------------------------------------------------*
#  def SUBPAR                                                  *:
#  Purpose: Dry matter partitioning fractions: leaves, stem and storage organs.
#----------------------------------------------------------------------*

# DEFINE_CALL  SUBPAR (INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                      OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

def SUBPAR (NPART, TRANRF, NNI, FRTWET, FLVT, FSTT, FSOT, FSHMOD):
  
  
  
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

    return FSHMOD, FRT, FLV, FST, FSO # FLVMOD removed from signature - WdW



# ---------------------------------------------------------------------*
#  def GROWTH                                                   *:
#  Purpose: To compute the total growth rate.                          *
# ---------------------------------------------------------------------*

# DEFINE_CALL GROWTH(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                     INPUT,                         OUTPUT,OUTPUT)




def GROWTH(TIME, DOYEM, DTR, K, NLUE, LAI, LUE, TRANRF, NNI):
    
    PARINT = 0.5 * DTR * (1.- exp(-K*LAI)) * INSW(TIME-DOYEM,0.,1.)
    
    if(TRANRF  <=  NNI):
        #  Water stress is more severe as compared to nitrogen stress and
        #  partitioning will follow the original assumptions of LINTUL2*
        
        GTOTAL = LUE * PARINT * TRANRF
    
    else:
    
        #  Nitrogen stress is more severe as compared to water stress and the
        #  less partitioning to leaves will go to the roots*
        
        GTOTAL = LUE * PARINT * exp(-NLUE*(1.0-NNI))
  
#     return PARINT, GTOTAL
    return GTOTAL



# ---------------------------------------------------------------------*
#  def RELGR                                                    *:
#  Purpose: To compute the relative growth rate of roots, leaves, stem *
#        and storage organs.                                        *
# ---------------------------------------------------------------------*

# DEFINE_CALL RELGR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                   INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                                         OUTPUT,OUTPUT,OUTPUT,OUTPUT)

def RELGR(TIME, DOYEM, EMERG, WLVGI, WRTLI, WSTI, WSOI, GTOTAL, 
          FLV, FRT, FST, FSO, DLV, DRRT, DELT):
  
  
    if (TIME  >= DOYEM  and  EMERG  ==  1.0):
        RWLVG = GTOTAL * FLV - DLV
        RWRT  = GTOTAL * FRT - DRRT
        RWST  = GTOTAL * FST
        RWSO  = GTOTAL * FSO
        
    else:
        RWLVG = 0.
        RWRT  = 0.
        RWST  = 0.
        RWSO  = 0.
  
    return RWLVG, RWRT, RWST, RWSO



# ---------------------------------------------------------------------*
#  def RNUSUB                                                   *:
#  Purpose: To compute the partitioning of the total N uptake rate     *
#        (NUPTR) over the leaves, stem, and roots.                  *
# ---------------------------------------------------------------------*

# DEFINE_CALL RNUSUB(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    OUTPUT,OUTPUT,OUTPUT)

def RNUSUB(TIME, DOYEM, EMERG, NDEML, NDEMS, NDEMR, NUPTR,         
           NDEMTO, ANLVI, ANSTI, ANRTI, DELT):
  
  
    if (TIME  >= DOYEM  and  EMERG  == 1):
        RNULV = (NDEML / notNull(NDEMTO))* NUPTR
        RNUST = (NDEMS / notNull(NDEMTO))* NUPTR
        RNURT = (NDEMR / notNull(NDEMTO))* NUPTR
    
    else:
        RNULV = 0.
        RNUST = 0.
        RNURT = 0.
    
    return RNULV, RNUST, RNURT



# ---------------------------------------------------------------------*
#  def NOPTM                                                    *:
#  Purpose: To compute the maximum nitrogen concentration of the organs*
# ---------------------------------------------------------------------*

# DEFINE_CALL NOPTM(INPUT,INPUT,INPUT,    OUTPUT,OUTPUT)

def NOPTM(FRNX, NMAXLV, NMAXST):
    
    NOPTLV= FRNX * NMAXLV
    NOPTST= FRNX * NMAXST
  
    return NOPTLV, NOPTST



# ---------------------------------------------------------------------*
#  def NDEMND                                                   *:
#  Purpose: To compute the nitrogen demand of the vegetative organs.   *
# ---------------------------------------------------------------------*

# DEFINE_CALL NDEMND(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
#                    INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    INPUT,INPUT,INPUT,INPUT,...
#                    OUTPUT,OUTPUT,OUTPUT,OUTPUT)

def NDEMND(NMAXLV, NMAXST, NMAXRT, NMAXSO, WLVG, WST, WRT, WSO,     
           RWLVG, RWST, RWRT, RWSO,                             
           ANLV, ANST, ANRT, ANSO, TCNT, DELT):
    
    NDEML  =  max (NMAXLV*WLVG -  ANLV, 0.)
    NDEMS  =  max (NMAXST*WST  - ANST, 0.)
    NDEMR  =  max (NMAXRT*WRT  - ANRT, 0.)
    NDEMSO =  max (NMAXSO*WSO  - ANSO, 0.)/TCNT
  
    return NDEML, NDEMS, NDEMR, NDEMSO




# ---------------------------------------------------------------------*
#  def NTRLOC                                                   *:
#  Purpose: To compute the translocatable N in the organs.             *
# ---------------------------------------------------------------------*

# DEFINE_CALL NTRLOC(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    INPUT,INPUT,         OUTPUT,OUTPUT,OUTPUT,OUTPUT)

def NTRLOC(ANLV, ANST, ANRT, WLVG, WST, WRT, RNFLV, RNFST, RNFRT, FNTRT):
  
    ATNLV = max (0. , ANLV-WLVG*RNFLV)
    ATNST = max (0. , ANST-WST*RNFST)
    ATNRT = min((ATNLV + ATNST) * FNTRT, ANRT-WRT*RNFRT)
    ATN   = ATNLV +  ATNST + ATNRT
    
    return ATNLV, ATNST, ATNRT, ATN



# ---------------------------------------------------------------------*
#  def NTRANS                                                   *:
#  Purpose: To compute the N translocated from different organs.       *
# ---------------------------------------------------------------------*

# DEFINE_CALL NTRANS(INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT,OUTPUT,OUTPUT)

def NTRANS(RNSO, ATNLV, ATNST, ATNRT, ATN):
    
    RNTLV= RNSO* ATNLV/ notNull(ATN)
    RNTST= RNSO* ATNST/ notNull(ATN)
    RNTRT= RNSO* ATNRT/ notNull(ATN)

    return RNTLV, RNTST, RNTRT



# ---------------------------------------------------------------------*
#  def RNLD                                                     *:
#  Purpose: To compute the N loss due to death of leaves and roots.    *
# ---------------------------------------------------------------------*

# DEFINE_CALL RNLD(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                                          OUTPUT,OUTPUT,OUTPUT)

def RNLD(DVS, WRT, RDRRT, RNFLV, DLV, RNFRT, DVSDR):
    
    DRRT= 0. if (DVS < DVSDR) else WRT * RDRRT
    
    RNLDLV= RNFLV* DLV
    RNLDRT= RNFRT* DRRT
  
    return DRRT, RNLDLV, RNLDRT



# ---------------------------------------------------------------------*
#  def NNINDX                                                   *:
#  Purpose: To compute Nitrogen Nutrition Index                        *
# ---------------------------------------------------------------------*

def NNINDX(TIME, DOYEM, EMERG, NFGMR, NRMR, NOPTMR):
  
    TINY=0.001
    
    if (TIME  >= DOYEM  and  EMERG  ==  1.):
        NNI= limit (TINY, 1.0, ((NFGMR-NRMR)/notNull(NOPTMR-NRMR)))
    else:
        NNI = 0.0
  
    return NNI



# ---------------------------------------------------------------------*
#  def DEATHL                                                   *:
#  Purpose: To compute the relative death rate of leaves due to age, *
#        shading amd due to nitrogen stress.                        *
# ---------------------------------------------------------------------*

# DEFINE_CALL DEATHL(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
#                    INPUT,INPUT,INPUT,INPUT,          OUTPUT,OUTPUT,...
#                       OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT,OUTPUT)

def DEATHL(TIME, DOYEM, TSUM, TSUMAG, RDRTMP, RDRSHM, LAI,         
           LAICR, WLVG, RDRNS, NNI, SLA):
  
    RDRDV = 0. if (TSUM  <  TSUMAG) else RDRTMP
    
    RDRSH = max(0., RDRSHM * (LAI-LAICR) / LAICR)
    RDR   = max(RDRDV, RDRSH)
    
    if (TIME  >= DOYEM  and  NNI  <  1.):
        DLVNS   = WLVG *RDRNS * (1.-NNI)
        DLAINS  = DLVNS * SLA
    else:
        DLVNS   = 0.
        DLAINS  = 0.
    
    DLVS  = WLVG * RDR
    DLAIS = LAI *  RDR
    
    DLV   = DLVS + DLVNS
    DLAI  = DLAIS + DLAINS

#     return RDRDV, RDRSH, RDR, DLV, DLVS, DLVNS, DLAIS, DLAINS, DLAI
    return DLV, DLAI

