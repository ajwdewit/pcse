"""
Copyright (c) 2004-2012 Alterra, Wageningen-UR

This module provides an implementation of the WOFOST crop model in the
python interpreter.

See Also
--------
* http://www.wageningenur.nl/wofost
* http://wofost.wikispaces.com
"""
__author__ = "Allard de Wit <allard.dewit@wur.nl>"
__license__ = "European Union Public License"
__version__ = "5.0.0"

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




