# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
""" Collection of tests for PCSE.
"""
import unittest
import warnings

from . import test_assimilation
from . import test_abioticdamage
from . import test_partitioning
from . import test_evapotranspiration
from . import test_respiration
from . import test_wofost72
from . import test_penmanmonteith
from . import test_agromanager
from . import test_wofost80
from . import test_lintul3

def make_test_suite(dsn=None):
    """Assemble test suite and return it
    """
    allsuites = unittest.TestSuite([test_abioticdamage.suite(),
                                    # test_assimilation.suite(),  # skip test because test inputs do not include TMIN
                                    test_partitioning.suite(),
                                    test_evapotranspiration.suite(),
                                    test_respiration.suite(),
                                    test_penmanmonteith.suite(),
                                    test_agromanager.suite(),
                                    test_wofost72.suite(dsn),
                                    # test_lintul3.suite(),
                                    test_wofost80.suite()
                                    ])
    return allsuites

def test_all(dsn=None):
    """Assemble test suite and run the test using the TextTestRunner
    """
    allsuites = make_test_suite(dsn)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        unittest.TextTestRunner(verbosity=2).run(allsuites)
