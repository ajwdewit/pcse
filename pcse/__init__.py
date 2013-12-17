"""
Copyright (c) 2004-2012 Alterra, Wageningen-UR

This module provides an implementation of the WOFOST crop model in the
python interpreter.

See Also
--------
* http://www.wofost.wur.nl
* http://wofost.wikispaces.com
"""
__author__ = "Allard de Wit <allard.dewit@wur.nl>"
__license__ = "European Union Public License"
__version__ = "5.0.0"

from . import examples
from . import util
from . import db
from . import fileinput
from . import tests
from . import agromanagement
from . import soil
from . import crop
from .pywofost import PyWofost
from .start_wofost import start_wofost

def test(dsn=None):
    """Run all available tests for PyWOFOST."""
    tests.test_all(dsn)
