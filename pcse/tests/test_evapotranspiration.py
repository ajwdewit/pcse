# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import unittest
from datetime import date

from ..crop.evapotranspiration import Evapotranspiration
from ..base_classes import VariableKiosk, ParameterProvider
from .test_data import pot_evtra_testdata, wl_evtra_testdata1,\
    wl_parvalue_dict1, wl_evtra_testdata2, wl_parvalue_dict2

class Test_PotentialEvapotranspiration(unittest.TestCase):
    
    test_vars = ["EVWMX","EVSMX","TRAMX","TRA"]
    def setUp(self):
        self.test_data = pot_evtra_testdata
        # settings for kiosk
        self.kiosk = VariableKiosk()
        self.kiosk.register_variable(0, "DVS", type="S", publish=True)
        self.kiosk.register_variable(0, "LAI", type="S", publish=True)
        self.kiosk.register_variable(0, "SM", type="S", publish=True)
        # Initialize ET module
        cropdata = {"CFET":1.00, "KDIFTB":[0., 0.6, 2.0, 0.6], "DEPNR":4.5,
                    "IOX":0, "IAIRDU":0, "CRAIRC":-99.}
        soildata = {"SM0":0.4, "SMFCF":0.3, "SMW":0.1}
        sitedata = timerdata = {}
        parvalues = ParameterProvider(sitedata, timerdata, soildata, cropdata)
        dummyday = date(2000,1,1)
        self.evtra = Evapotranspiration(dummyday, self.kiosk, parvalues)
    
    def runTest(self):
        day = date(2000,1,1) # dummy date value
        SM = 0.3
        for drvref in self.test_data:
            self.kiosk.set_variable(0, "LAI", drvref.LAI)
            self.kiosk.set_variable(0, "DVS", drvref.DVS)
            self.kiosk.set_variable(0, "SM", SM)
            self.evtra(day, drvref)
            for testvar in self.test_vars:
                test_value = self.evtra.get_variable(testvar)
                ref_value  = getattr(drvref, testvar)
                self.assertTrue(abs(test_value - ref_value) < 0.0001)


class Test_WaterLimitedEvapotranspiration1(unittest.TestCase):
    """Testing Wat-lim evapotranspiration with DEPNR=4.5
    """
    
    test_vars = ["EVWMX","EVSMX","TRAMX","TRA"]
    def setUp(self):
        self.test_data = wl_evtra_testdata1
        # settings for kiosk
        self.kiosk = VariableKiosk()
        self.kiosk.register_variable(0, "DVS", type="S", publish=True)
        self.kiosk.register_variable(0, "LAI", type="S", publish=True)
        self.kiosk.register_variable(0, "SM", type="S", publish=True)
        # Initialize ET module
        cropdata = {}
        for key in ["CFET", "KDIFTB", "DEPNR","IOX","IAIRDU", "CRAIRC"]:
            cropdata[key] = wl_parvalue_dict1[key]
        soildata = {}
        for key in ["SM0", "SMFCF", "SMW"]:
            soildata[key] = wl_parvalue_dict1[key]
        sitedata = timerdata = {}
        parvalues = ParameterProvider(sitedata, timerdata, soildata, cropdata)
        dummyday = date(2000, 1, 1)
        self.evtra = Evapotranspiration(dummyday, self.kiosk, parvalues)
    
    def runTest(self):
        day = date(2000,1,1) # dummy date value
        for drvref in self.test_data:
            self.kiosk.set_variable(0, "LAI", drvref.LAI)
            self.kiosk.set_variable(0, "DVS", drvref.DVS)
            self.kiosk.set_variable(0, "SM", drvref.SM)
            self.evtra(day, drvref)
            for testvar in self.test_vars:
                test_value = self.evtra.get_variable(testvar)
                ref_value  = getattr(drvref, testvar)
                self.assertTrue(abs(test_value - ref_value) < 0.0001)


class Test_WaterLimitedEvapotranspiration2(unittest.TestCase):
    """Testing Wat-lim evapotranspiration with DEPNR=2.5
    """
    
    test_vars = ["EVWMX","EVSMX","TRAMX","TRA"]
    def setUp(self):
        self.test_data = wl_evtra_testdata2
        # settings for kiosk
        self.kiosk = VariableKiosk()
        self.kiosk.register_variable(0, "DVS", type="S", publish=True)
        self.kiosk.register_variable(0, "LAI", type="S", publish=True)
        self.kiosk.register_variable(0, "SM", type="S", publish=True)
        # Initialize ET module
        cropdata = {}
        for key in ["CFET", "KDIFTB", "DEPNR","IOX","IAIRDU", "CRAIRC"]:
            cropdata[key] = wl_parvalue_dict2[key]
        soildata = {}
        for key in ["SM0", "SMFCF", "SMW"]:
            soildata[key] = wl_parvalue_dict2[key]
        sitedata = timerdata = {}
        parvalues = ParameterProvider(sitedata, timerdata, soildata, cropdata)
        dummyday = date(2000,1,1)
        self.evtra = Evapotranspiration(dummyday, self.kiosk, parvalues)
    
    def runTest(self):
        day = date(2000,1,1) # dummy date value
        for drvref in self.test_data:
            self.kiosk.set_variable(0, "LAI", drvref.LAI)
            self.kiosk.set_variable(0, "DVS", drvref.DVS)
            self.kiosk.set_variable(0, "SM", drvref.SM)
            self.evtra(day, drvref)
            for testvar in self.test_vars:
                test_value = self.evtra.get_variable(testvar)
                ref_value  = getattr(drvref, testvar)
                self.assertTrue(abs(test_value - ref_value) < 0.0001)

            
def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_PotentialEvapotranspiration))
    suite.addTest(unittest.makeSuite(Test_WaterLimitedEvapotranspiration1))
    suite.addTest(unittest.makeSuite(Test_WaterLimitedEvapotranspiration2))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())

