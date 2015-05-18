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
 
# Model parameters
# ----------------
DOYEM  = 90.
DRATE  = 30.
DVSDR  = 1.0
DVSNLT = 1.0
DVSNT  = 0.8
FNTRT  = 0.15
FRNX   = 0.5
IRRIGF = 1.0
K      = 0.6
LAICR  = 4.0
LRNR   = 0.50
LSNR   = 0.50
LUE    = 2.8
NFRLVI = 0.06
NFRRTI = 0.03
NFRSTI = 0.03
NLAI   = 1.0
NLUE   = 0.20
NMAXSO = 0.0165
NPART  = 1.0
NSLA   = 1.0
RDRSHM = 0.03
RGRL   = 0.009
RNFLV  = 0.004
RNFRT  = 0.002
RNFST  = 0.002
ROOTDM = 1.2
RRDMAX = 0.012
SLAC   = 0.022
TBASE  = 0.		# TBASEM = Float(-99.)  # Base temp. for emergence
TCNT   = 10.
TRANCO = 8.
TSUMAG = 800.
TSUMAN = 800.	# TSUM1  = Float(-99.)  # Temperature sum emergence to anthesis
TSUMMT = 1030.	# TSUM2  = Float(-99.)  # Temperature sum anthesis to maturity
WCAD   = 0.10
WCFC   = 0.40
WCI    = 0.40
WCST   = 0.50
WCSUBS = 0.30
WCWET  = 0.45
WCWP   = 0.20
WMFAC  = 0.0
	
# Interpolation functions
# -----------------------
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

SLACF = [
    0.0, 1.0,
    2.0, 1.0,
    2.1, 1.0 ]

NMXLV = [
    0.0, 0.06,
    0.4, 0.04,
    0.7, 0.03,
    1.0, 0.02,
    2.0, 0.014,
    2.1, 0.014 ]

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

FERTAB = [
    0.0, 0.0,
    99., 0.0,
    100., 10.0,
    101., 0.0,
    124.0, 0.0,
    125., 5.0,
    126., 0.0,
    200., 0.0 ]

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
SCALExxx = [0.]*27

# lines below inserted to have at least some realistic parameters for Phenology - WdW	
# -----------------------------------------------------------------------------------
# ** Rice is transplanted, no simulation before emergence.
TBASEM = 00.0 	# lower threshold temp. for emergence [cel]
TEFFMX = 00.0 	# max. eff. temp. for emergence [cel]
TSUMEM = 00. 	# temperature sum from sowing to emergence [cel d]
# ** phenology
IDSL = 0 		# indicates whether pre-anthesis development depends
				# on temp. (=0), daylength (=1) , or both (=2)
DLO = 11.5 		# optimum daylength for development [hr]
DLC = 0.0 		# critical daylength (lower threshold) [hr]
TSUM1 =1420. 	# temperature sum from emergence to anthesis [cel d]
TSUM2 = 680. 	# temperature sum from anthesis to maturity [cel d]
DTSMTB = [0.00, 0.00, 	# daily increase in temp. sum
          8.00, 0.00, 	# as function of av. temp. [cel; cel d]
          30.00, 24.00,
          42.50, 0.00]
DVSI = 0.187 	# development stage start simulation (after transplanting)
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
 
