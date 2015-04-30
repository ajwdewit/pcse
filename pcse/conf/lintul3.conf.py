# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""PCSE configuration file for testing lintul3.

This configuration file defines the crop.
"""

from pcse.lintul.lintul3 import Lintul3
from pcse.empty_model import EmptyModel
from pcse.agromanagement import AgroManagementSingleCrop

# Module to be used for water balance
SOIL = EmptyModel

# Module to be used for the crop simulation itself
CROP = Lintul3

# Module to use for AgroManagement actions
AGROMANAGEMENT = AgroManagementSingleCrop

# variables to save at OUTPUT signals
# Set to an empty list if you do not want any OUTPUT
OUTPUT_VARS = ["TIME", "WAI", "DVS", "TSUM", "TAGBM", "WST", "WLVG", "WLVD", "WSO", "LAI", "NTAC", "WRT", 
               "GTSUM", "CBALAN", "TRANRF", "NNI", "SLA", "FRACT", "FRTWET", "FLVT", "FSTT", "FSOT", 
               "rWLVG", "rWST", "rWRT", "rWSO", "CUMPAR", "LUECAL", "NUPTT", "TTRAN", "TEVAP", "PEVAP", 
               "NBALAN", "WATBAL", "NUPTR", "TNSOIL", "NDEMTO", "RNSOIL", "FERTN", "FERTNS", "WA", 
               "TIRRIG", "TRAIN", "TEXPLO", "TRUNOF", "TDRAIN"]


# interval for OUTPUT signals, either "daily"|"dekadal"|"monthly"                                    
# For daily output you change the number of days between successive
# outputs using OUTPUT_INTERVAL_DAYS. For dekadal and monthly
# output this is ignored.
OUTPUT_INTERVAL = "daily"
OUTPUT_INTERVAL_DAYS = 1
OUTPUT_WEEKDAY = 0

# variables to save at SUMMARY_OUTPUT signals
# Set to an empty list if you do not want any SUMMARY_OUTPUT
SUMMARY_OUTPUT_VARS = []
