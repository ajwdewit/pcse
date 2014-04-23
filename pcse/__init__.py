# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""
Copyright (c) 2004-2014 Alterra, Wageningen-UR

The Python Crop Simulation Environment (PCSE) has been developed
to facilitate implementing crop simulation models that were 
developed in Wageningen. PCSE provides a set of building blocks
that on the one hand facilitate implementing the crop simulation 
models themselves and other hand allow to interface these models with 
external inputs and outputs (files, databases, webservers) 

PCSE builds on existing ideas implemented in the FORTRAN
Translator (FST) and its user interface FSE. It inherits ideas
regarding the rigid distinction between rate calculation
and state integration and the initialization of parameters
in a PCSE model. Moreover PCSE provides support for reusing
input files and weather files that are used by FST models.

PCSE currently provides an implementation of the WOFOST crop 
simulation model and variants of WOFOST with extended
capabilities.

See Also
--------
* http://www.wageningenur.nl/wofost
* http://wofost.wikispaces.com
"""
__author__ = "Allard de Wit <allard.dewit@wur.nl>"
__license__ = "European Union Public License"
__version__ = "0.9.0"

import sys, os

# First define and run setup before importing the rest of the stuff
def setup():
    """
    Set up the .pcse folder and user settings file, add ~/.pcse to the
    sys.path and setup logging.
    """

    user_home = os.path.expanduser("~")
    pcse_user_home = os.path.join(user_home, ".pcse")
    if not os.path.exists(pcse_user_home):
        os.mkdir(pcse_user_home)

    # Add PCSE home to python PATH
    sys.path.append(pcse_user_home)

    # Check existence of user settings file. If not exists, create it.
    user_settings_file = os.path.join(pcse_user_home, "user_settings.py")
    if not os.path.exists(user_settings_file):
        pcse_dir = os.path.dirname(__file__)
        default_settings_file = os.path.join(pcse_dir, "settings", "default_settings.py")
        lines = open(default_settings_file).readlines()
        with open(user_settings_file, "w") as fp:
            for line in lines:
                if line.startswith(("#", '"', "'", "import")):
                    cline = line
                elif len(line.strip()) == 0: # empty line
                    cline = line
                else:
                    cline = "# " + line
                fp.write(cline)

setup()

import logging.config
from .settings import settings
logging.config.dictConfig(settings.LOG_CONFIG)

from . import examples
from . import util
from . import db
from . import fileinput
from . import tests
from . import agromanagement
from . import soil
from . import crop
from .start_wofost import start_wofost

def test(dsn=None):
    """Run all available tests for PCSE."""
    tests.test_all(dsn)




