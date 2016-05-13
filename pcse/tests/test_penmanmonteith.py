# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import os
import unittest
import csv
import math
from datetime import datetime

from ..util import penman_monteith


# Lat, lon, AngstA, AngstB for each station
stations = {8011: (43.5667, -6.0442, 130., 0.243, 0.541),
            2297: (64.4833, 21.5833, 37., 0.241, 0.419),
            10870:(48.3667, 11.8, 453., 0.225, 0.505),
            8443: (36.75, -5.1667, 3., 0.3064, 0.4846)}

# Allow an ET0 difference of 0.2 mm to account for limited precision
# in CGMS database.
delta = 0.2

class DottedContainer(object):
    def __init__(self, d):
        for key, value in d.items():
            if key == 'DAY':
                value = (datetime.strptime(value, '%d-%b-%y').date())
            elif key == 'STATION_NUMBER':
                value = int(value)
            elif key == 'WINDSPEED_10M':
                z = 10.
                u10 = float(value)
                value = (u10 * 4.87) / math.log((67.8*z) - 5.42)
                key = "WINDSPEED_2M"
            elif key == 'RADIATION':
                value = float(value)*1000.
            else:
                value = float(value)
            setattr(self, key, value)

class PenmanMonteith_TestingTemplate(unittest.TestCase):
    test_input = None # Must be replaced by real csv filename in subclass

    def setUp(self):
        test_dir = os.path.dirname(__file__)
        csv_file = os.path.join(test_dir, "test_data", self.test_input)
        self.csv_data = []
        with open(csv_file) as fp:
            dct = csv.Sniffer().sniff(fp.read(1024))
            fp.seek(0)
            reader = csv.DictReader(fp, dialect=dct)
            for row in reader:
                try:
                    self.csv_data.append(DottedContainer(row))
                except (TypeError, ValueError):
                    continue

    def runTest(self):
        LAT = None
        for rec in self.csv_data:
            if LAT is None:
                LAT, LON, ELEV, ANGSTA, ANGSTB = stations[rec.STATION_NUMBER]
            ET0 = penman_monteith(rec.DAY, LAT, ELEV, rec.TEMPERATURE_MIN,
                                  rec.TEMPERATURE_MAX, rec.RADIATION,
                                  rec.VAPOURPRESSURE, rec.WINDSPEED_2M)
            self.assertLess(abs(ET0 - rec.ET0), delta)

class Test_PenmanMonteith1(PenmanMonteith_TestingTemplate):
    test_input = '8011_asturias.csv'

class Test_PenmanMonteith2(PenmanMonteith_TestingTemplate):
    test_input = '2297_sweden.csv'

class Test_PenmanMonteith3(PenmanMonteith_TestingTemplate):
    test_input = '8443_ronda.csv'

class Test_PenmanMonteith4(PenmanMonteith_TestingTemplate):
    test_input = '10870_munchen_flughafen.csv'

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_PenmanMonteith1))
    suite.addTest(unittest.makeSuite(Test_PenmanMonteith2))
    suite.addTest(unittest.makeSuite(Test_PenmanMonteith3))
    suite.addTest(unittest.makeSuite(Test_PenmanMonteith4))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())
