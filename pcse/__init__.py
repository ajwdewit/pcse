# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""
Copyright (c) 2004-2014 Alterra, Wageningen-UR

The Python Crop Simulation Environment (PCSE) has been developed
to facilitate implementing crop simulation models that were 
developed in Wageningen. PCSE provides a set of building blocks
that on the one hand facilitates implementing the crop simulation
models themselves and other hand allows to interface these models with
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
from __future__ import print_function
__author__ = "Allard de Wit <allard.dewit@wur.nl>"
__license__ = "European Union Public License"
__version__ = "5.0.0"
__stable__ = False

import sys, os

# First define and run setup before importing the rest of the stuff
def setup():
    """
    Set up the .pcse folder and user settings file, add ~/.pcse to the
    sys.path.
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

# from . import examples
from . import util
from . import db
from . import fileinput
from . import tests
from . import agromanagement
from . import soil
from . import crop
from .start_wofost import start_wofost

# If no PCSE demo database, build it!
pcse_db_file = os.path.join(settings.PCSE_USER_HOME, "pcse.db")
if not os.path.exists(pcse_db_file):
    print("Building PCSE demo database at: %s ..." % pcse_db_file, end=" ")
    pcse_home = os.path.dirname(__file__)
    pcse_db_dump_file = os.path.join(pcse_home, "db", "pcse", "pcse_dump.sql")
    try:
        util.load_SQLite_dump_file(pcse_db_dump_file, pcse_db_file)
        print("OK")
    except Exception as e:
        logger = logging.getLogger()
        msg1 = "Failed to create the PCSE demo database: %s" % e
        msg2 = "PCSE will likely be functional, but some tests and demos may fail."
        logger.warn(msg1)
        logger.warn(msg2)

if not __stable__:
    print("Warning: You are running a PCSE development version:  %s" % __version__)

def test(dsn=None):
    """Run all available tests for PCSE."""
    tests.test_all(dsn)




