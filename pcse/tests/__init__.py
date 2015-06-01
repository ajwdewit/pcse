# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
""" Collection of tests for PyWOFOST and individual components.
"""
import unittest
from . import test_assimilation
from . import test_abioticdamage
from . import test_partitioning
from . import test_evapotranspiration
from . import test_respiration
from . import test_wofost
from . import test_penmanmonteith

def make_test_suite(dsn=None):
    """Assemble test suite and return it
    """
    allsuites = unittest.TestSuite([test_abioticdamage.suite(),
                                    test_assimilation.suite(),
                                    test_partitioning.suite(),
                                    test_evapotranspiration.suite(),
                                    test_respiration.suite(),
                                    test_penmanmonteith.suite(),
                                    test_wofost.suite(dsn)])
    return allsuites

def test_all(dsn=None):
    """Assemble test suite and run the test using the TextTestRunner
    """
    allsuites = make_test_suite(dsn)
    unittest.TextTestRunner(verbosity=2).run(allsuites)
