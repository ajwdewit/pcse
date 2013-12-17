""" Collection of tests for PyWOFOST and individual components.
"""
import unittest
import test_assimilation
import test_abioticdamage
import test_partitioning
import test_evapotranspiration
import test_respiration
import test_pywofost

def test_all(dsn=None):
    allsuites = unittest.TestSuite([test_abioticdamage.suite(), 
                                    test_assimilation.suite(),
                                    test_partitioning.suite(),
                                    test_evapotranspiration.suite(),
                                    test_respiration.suite(),
                                    test_pywofost.suite(dsn)])
    unittest.TextTestRunner(verbosity=2).run(allsuites)
