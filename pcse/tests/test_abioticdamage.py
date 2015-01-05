# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import unittest
from datetime import date

from ..crop.abioticdamage import FROSTOL
from ..base_classes import VariableKiosk
from .test_data import frostol_testdata

#------------------------------------------------------------------------------
class Test_FROSTOL(unittest.TestCase):
    """Testing suite for FROSTOL.
    
    Test data kindly provided by Anne Kari Bergjord.
    """

    test_vars = ["LT50T","RH","RDH_TEMP","RDH_RESP","RDH_TSTR"]
    def setUp(self):
        self.testdata = frostol_testdata
        # get parameter values from first testdata record
        r = self.testdata[1]
        parvalues = {"LT50C":r.LT50C,
                    "IDSL":2,
                    "FROSTOL_D":r.FROSTOL_D,
                    "FROSTOL_H":r.FROSTOL_H,
                    "FROSTOL_R":r.FROSTOL_R,
                    "FROSTOL_S":r.FROSTOL_S,
                    "FROSTOL_SDBASE":0.,
                    "FROSTOL_SDMAX":12.5,
                    "FROSTOL_KILLCF":1.019,
                    "ISNOWSRC":1,
                    "CROWNTMPA":0.5,
                    "CROWNTMPB":0.2}
        # Setup variable kiosk and register variables
        self.kiosk = VariableKiosk()
        self.kiosk.register_variable(0, "ISVERNALISED", type="S", publish=True)
        self.kiosk.register_variable(0, "SNOWDEPTH", type="S", publish=True)
        # Initialize FROSTOL
        dummyday = date(2000,1,1)
        self.frostol = FROSTOL(dummyday, self.kiosk, parvalues, testing=True)

    #@unittest.skip("FROSTOL test failing because of problem with test")
    def runTest(self):
        for day in range(1, 252):
            # reference data and driving variables
            drvref = self.testdata[day]

            # Set values in kiosk
            vern = False if (drvref.fV < 0.99) else True
            self.kiosk.set_variable(0, "ISVERNALISED", vern)
            self.kiosk.set_variable(0, "SNOWDEPTH", drvref.snow_depth)

            # calculated rates
            self.frostol.calc_rates(day, drvref)
            
            # Assert simulation results almost equal to reference values
            for var in self.test_vars:
                refvalue = getattr(drvref, var)
                modvalue = self.frostol.get_variable(var)
                self.assertTrue(abs(refvalue - modvalue) < 0.5)

            # integrate
            self.frostol.integrate(day)

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_FROSTOL))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())

