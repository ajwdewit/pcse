# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

import unittest
from datetime import date

from ..base_classes import VariableKiosk
from ..crop.respiration import WOFOST_Maintenance_Respiration
from .test_data import respiration_testdata


class Test_WOFOSTMaintenanceRespiration(unittest.TestCase):
    """Unit tests for WOFOST maintenance respiration"""
    def setUp(self):
        self.testdata = respiration_testdata
        self.cropdata = {"RMR":0.015,"RML":0.03,"RMS":0.015,"RMO":0.01,
                        "RFSETB":[0.0,1.0,2.0,1.0],"Q10":2.0}
        self.kiosk = VariableKiosk()
        self.kiosk.register_variable(0, "WRT", type="S", publish=True)
        self.kiosk.register_variable(0, "WLV", type="S", publish=True)
        self.kiosk.register_variable(0, "WST", type="S", publish=True)
        self.kiosk.register_variable(0, "WSO", type="S", publish=True)
        self.kiosk.register_variable(0, "DVS", type="S", publish=True)
        # initialize respiration model
        dummyday = date(2000,1,1)
        self.respiration = WOFOST_Maintenance_Respiration(dummyday, self.kiosk,
                                                        self.cropdata)

    def runTest(self):
        day = date(2000,1,1) # dummy date value
        for drvref in self.testdata:
            self.kiosk.set_variable(0, "WRT", drvref.WRT)
            self.kiosk.set_variable(0, "WLV", drvref.WLV)
            self.kiosk.set_variable(0, "WST", drvref.WST)
            self.kiosk.set_variable(0, "WSO", drvref.WSO)
            self.kiosk.set_variable(0, "DVS", drvref.DVS)
            PMRES = self.respiration(day, drvref)
            self.assertTrue(abs(PMRES - drvref.PMRES) < 0.0001)

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_WOFOSTMaintenanceRespiration))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())


