#----------------------------------------------------------------------#
# Information about this file                                          #
# ===========================                                          #
# Contents      : Model data file                                      #
# Creator       : FST translator version 4.00                          #
# Creation date : 16-Apr-2015, 11:44:49                                #
# Source file   : LINTUL3-SPRINGWHEAT-TEST.FST                         #
#----------------------------------------------------------------------#
 
# contains:
# - Initial constants as far as specified with INCON statements,
# - Model parameters,
# - AFGEN/CSPLIN interpolation functions,
# - A SCALE array in case of a general translation
 
# Initial constants
# -----------------
ZERO   = 0.
TSUMI  = 0.
WLVGI  = 2.4
WSTI   = 0.0
WRTLI  = 3.6
WSOI   = 0.0
ROOTDI = 0.10

# Soil hydraulic properties
WCAD   = 0.10 			
WCWP   = 0.20
WCFC   = 0.40
WCWET  = 0.45
WCST   = 0.50

DOYEM  = 90.  			# The day of the year on which crop emerges.
DRATE  = 30.  			# Maximum drainage rate of the soil (mm/day) 
DVSDR  = 1.0  			# Development stage above which death of leaves and roots start.
DVSNLT = 1.0  			# development stage N-limit
DVSNT  = 0.8  			# development stage N-threshold
FNTRT  = 0.15 			# Nitrogen translocation from roots as a fraction of the total amount of nitrogen translocated from leaves and stem.
FRNX   = 0.5  			# Optimal N concentration as the fraction of maximum N concentration.
K      = 0.6  			# light extinction coefficient
LAICR  = 4.0  			# (oC d)-1, critical LAI above which mutual shading of leaves occurs,
LRNR   = 0.50			
LSNR   = 0.50			
LUE    = 2.8  			# Light use efficiency.
NFRLVI = 0.06 			# Initial fraction of N (g N g-1 DM) in leaves.
NFRRTI = 0.03 			# Initial fraction of N (g N g-1 DM) in roots.
NFRSTI = 0.03 			# Initial fraction of N (g N g-1 DM) in stem.
NLAI   = 1.0  			# Coefficient for the effect of N stress on LAI reduction(during juvenile phase)
NLUE   = 0.20 			# Extinction coefficient for  Nitrogen distribution down the canopy
NMAXSO = 0.0165
NPART  = 1.0  			# Coefficient for the effect of N stress on leaf biomass reduction
NSLA   = 1.0  			# Coefficient for the effect of N stress on SLA reduction
RDRRT  = 0.03 			# Relative death rate of roots.
RDRSHM = 0.03 			# and the maximum relative death rate of leaves due to shading.
RGRL   = 0.009			# Relative growth rate of LAI at the exponential growth phase
RNFLV  = 0.004			# Residual N concentration in leaves
RNFRT  = 0.002			# Residual N concentration in roots.
RNFST  = 0.002			# Residual N concentration in stem
ROOTDM = 1.2  			# Maximum root depth for a rice crop.
RRDMAX = 0.012			# Maximum rate of increase in rooting depth (m d-1) for a rice crop.
SLAC   = 0.022			# Specific leaf area constant.
TBASE  = 0.   			# Base temperature for spring wheat crop.
TCNT   = 10.  			# Time coefficient(days) for N translocation.
TRANCO = 8.   			# Transpiration constant (mm/day) indicating the level of drought tolerance of the wheat crop.
TSUMAG = 800. 			# Temperature sum for ageing of leaves
TSUMAN = 800. 			# Temperature sum for anthesis  [corresponds to TSUM1 = Float(-99.)# Temperature sum emergence to anthesis]
TSUMMT = 1030.			# Temperature sum for maturity  [corresponds to TSUM2= Float(-99.) # Temperature sum anthesis to maturity
WCI    = 0.40 			# Initial water content in cm3 of water/(cm3 of soil).
WCSUBS = 0.30 			# water content subsoil (?)

# Water management
# We will represent two water management situations in the model
# as irrigated up to the field capacity:                  WMFAC= False
# as irrigated up to saturation, thus mimicking flooding; WMFAC= True
# In both cases parameter IRRIGF must be taken 1.0. "Irrigation" thus
WMFAC  = False
IRRIGF = True			# Irrigation factor 

# Interpolation functions
# -----------------------
# Function to include the effect of photoperiodicity
PHOTTB = [
    0.0, 0.0,
    8., 1.0,
    10., 1.0,
    12., 1.0,
    18., 1.0 ]

RDRT = [
    -10., 0.00,
    10., 0.02,
    15., .03,
    30., 0.05,
    50., 0.09 ]

# Leaf area correction function as a function of development stage, DVS.
SLACF = [
    0.0, 1.0,
    2.0, 1.0,
    2.1, 1.0 ]

#   Maximum N concentration in the leaves, from which the N-conc.values of the
#   stem and roots are derived, as a function of development stage.
NMXLV = [
    0.0, 0.06,
    0.4, 0.04,
    0.7, 0.03,
    1.0, 0.02,
    2.0, 0.014,
    2.1, 0.014 ]

# ********** Partitioning coefficients ***********************************
FRTTB = [
    0.0, 0.60,
    0.33, 0.58,
    0.40, 0.55,
    0.80, 0.10,
    1.00, 0.00,
    2.00, 0.0 ]

FLVTB = [
    0.0, 0.40,
    0.33, 0.42,
    0.40, 0.405,
    0.80, 0.36,
    1.00, 0.10,
    1.01, 0.00,
    2.00, 0.00 ]

FSTTB = [
    0.0, 0.00,
    0.33, 0.00,
    0.40, 0.045,
    0.80, 0.54,
    1.00, 0.90,
    1.01, 0.25,
    2.00, 0.00 ]

FSOTB = [
    0.0, 0.00,
    0.33, 0.00,
    0.40, 0.00,
    0.80, 0.00,
    1.00, 0.00,
    1.01, 0.75,
    2.00, 1.00 ]

# Fertilizer application as a function of TIME (g N m-2).
FERTAB = [
    0.0, 0.0,
    99., 0.0,
    100., 10.0,
    101., 0.0,
    124.0, 0.0,
    125., 5.0,
    126., 0.0,
    200., 0.0 ]

# Fertilizer nitrogen recovery fraction
NRFTAB = [
    0.0, 0.70,
    100., 0.70,
    125., 0.70,
    150., 0.70,
    200., 0.70 ]

 
# Size scale of state variables, used by the driver RKDRIV
# --------------------------------------------------------
# Zero's lead to relative integration errors EPS for states above 1.0
# and absolute integration errors EPS for states below 1.0. The value
# of EPS may be specified in a TIMER statement and is written to the
# file TIMER.DAT. The default value is 1.0E-4. The number of elements
# of SCALE is equal to the number of state variables.
# A non-zero value of SCALE causes Integration errors for the corresponding
# state variable to be evaluated relative to the value of SCALE.
SCALExxx = [0.]#27

# lines below inserted to have at least some realistic parameters for Phenology - WdW	
# -----------------------------------------------------------------------------------
# ## Rice is transplanted, no simulation before emergence.
TBASEM = 00.0 	# lower threshold temp. for emergence [cel]
TEFFMX = 00.0 	# max. eff. temp. for emergence [cel]
TSUMEM = 00. 	# temperature sum from sowing to emergence [cel d]
#  phenology
IDSL = 1 		# indicates whether pre-anthesis development depends
				# on temp. (=0), daylength (=1) , or both (=2)
DLO = 8. 		# optimum daylength for development [hr]
DLC = 0.0 		# critical daylength (lower threshold) [hr]
TSUM1 =TSUMAN 	# temperature sum from emergence to anthesis [cel d]
TSUM2 =TSUMMT 	# temperature sum from anthesis to maturity [cel d]
DTSMTB = [0.00, 0.00, 	# daily increase in temp. sum
          50.00, 50.00]	# as function of av. temp. [cel; cel d]
DVSI = 0 	# development stage start simulation (after transplanting)
DVSEND = 2.00 	# development stage at harvest or at
				# physiological maturity (= 2.0 at maturity [-])
# ------------- (end of insert) -----------------------------------------------------
 
# lines below inserted to have synonym parameters for Partitioning - WdW	
# -----------------------------------------------------------------------------------
FRTB     = FRTTB
FSTB     = FSTTB
FLTB     = FLVTB
FOTB     = FSOTB
# ------------- (end of insert) -----------------------------------------------------
 
#  Relative death rate of leaves due to N stress.
RDRNS  = 0.03

#  Relative death rate of roots.
RDRRT = 0.03

