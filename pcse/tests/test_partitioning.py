# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import unittest
from datetime import date

from ..crop.partitioning import DVS_Partitioning
from ..base_classes import VariableKiosk

#----------------------------------------------------------------------------
class Test_DVS_Partitioning(unittest.TestCase):
    """Unit test for DVS_Partitioning class.
    """
    DVS_ref = [0.1, 0.5, 1, 1.2, 1.5, 2, 2.5, -0.5]
    FR_ref = {-0.5: 0.5, 0.1: 0.5, 0.5: 0.12999999523162842,
              1: 0.019999999552965164, 1.2: 0.0, 1.5: 0.0, 2: 0.0, 2.5: 0.0}
    FL_ref = {0.5: 0.5, 1: 0.0, 2: 0.0, 0.1: 0.6499999761581421, 1.2: 0.0,
              -0.5: 0.6499999761581421, 2.5: 0.0, 1.5: 0.0}
    FS_ref = {0.5: 0.5, 1: 0.0, 2: 0.0, 0.1: 0.3499999940395355, 1.2: 0.0,
              -0.5: 0.3499999940395355, 2.5: 0.0, 1.5: 0.0}
    FO_ref = {0.5: 0.0, 1: 1.0, 2: 1.0, 0.1: 0.0, 1.2: 1.0, -0.5: 0.0,
              2.5: 1.0, 1.5: 1.0}

    def setUp(self):
        FLTB = [0,0.65, 0.1,0.65, 0.25,0.7, 0.5,0.5, 0.646,0.3,
                0.95,0 ,1,0, 2,0 ,0,0, 0,0]
        FOTB = [0,0, 0.1,0, 0.25,0, 0.5,0, 0.646,0, 0.95,0, 1,1,
                2,1, 0,0, 0,0]
        FRTB = [0,0.5, 0.1,0.5, 0.2,0.4, 0.35,0.22, 0.4,0.17, 0.5,0.13,
                0.7,0.07, 0.9,0.03, 1.2,0, 2,0]
        FSTB = [0,0.35, 0.1,0.35, 0.25,0.3, 0.5,0.5, 0.646,0.7,
                0.95,1, 1,0, 2,0, 0,0, 0,0]
        cropdata = {"FRTB":FRTB, "FLTB":FLTB, "FSTB":FSTB, "FOTB":FOTB}
        self.kiosk = VariableKiosk()
        # Register DVS with the kiosk
        self.kiosk.register_variable(1, "DVS", type="S", publish=True)
        self.kiosk.set_variable(1, "DVS", 0.)
        dummyday = date(2000,1,1)
        self.partitioning = DVS_Partitioning(dummyday, self.kiosk, cropdata)

    def runTest(self):
        for dvs in self.DVS_ref:
            self.kiosk.set_variable(1, "DVS", dvs)
            self.partitioning.integrate(None)
            (FR, FL, FS, FO) = self.partitioning.calc_rates(None, None)
            self.assertAlmostEqual(FR, self.FR_ref[dvs], 5)
            self.assertAlmostEqual(FL, self.FL_ref[dvs], 5)
            self.assertAlmostEqual(FS, self.FS_ref[dvs], 5)
            self.assertAlmostEqual(FO, self.FO_ref[dvs], 5)

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_DVS_Partitioning))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())
