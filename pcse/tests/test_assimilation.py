# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import unittest
from datetime import date

from ..base_classes import VariableKiosk
from ..crop.assimilation import WOFOST_Assimilation

#-------------------------------------------------------------------------------
class Test_WOFOST_Assimilation(unittest.TestCase):
    ref_results = [
        (10, 37.639999, 7.55583569E-02, 10598000., 4.9428568, 13.125000, 0.25510019, 36.513809),
        (20, 37.639999, 0.15837558, 13064000., 5.7571425, 13.150001, 0.45266339, 63.450977),
        (30, 37.639999, 0.24197641, 12981000., 5.1857142, 16.025000, 0.88996553, 123.85042),
        (40, 37.639999, 0.35700494, 15789000., 7.3000002, 18.525000, 1.8993067, 215.95074),
        (50, 37.639999, 0.48448902, 17403000., 8.0857143, 19.250000, 3.4342484, 301.07986),
        (60, 37.639999, 0.61263305, 19372000., 9.6999998, 20.350000, 4.9810600, 360.32208),
        (70, 37.639999, 0.75484121, 21821000., 10.800000, 24.775000, 5.9800878, 407.55038),
        (80, 37.639999, 0.90179235, 22034000., 9.8000002, 20.724998, 6.2283578, 431.77167),
        (90, 37.639999, 1.0259552, 16819000., 9.3000002, 15.974999, 5.9964075, 435.98254),
        (100, 37.639999, 1.1575972, 15669000., 10.157143, 18.825001, 5.8061633, 429.88681),
        (110, 37.639999, 1.2953160, 21771000., 12.214286, 21.400000, 4.5933480, 479.97336),
        (120, 37.639999, 1.4405867, 12674000., 10.128572, 17.775002, 2.8579102, 288.33484),
        (130, 37.639999, 1.5955797, 20809000., 14.014285, 19.075001, 1.2674762, 185.61636),
        (140, 37.639999, 1.7846010, 26069000., 16.585712, 27.300001, 0.29815927, 29.025328),
        (150, 37.639999, 1.9740561, 31392000., 16.357143, 29.799999, 0.0000000, 0.0000000)]
    
    class container(object):
        pass

    def setUp(self):
        cropdata = {"AMAXTB":[0,35.83, 1,35.83, 1.3,35.83, 2,4.48],
                    "EFFTB":[0.,0.45, 40.,0.45],
                    "KDIFTB":[0,0.6, 2,0.6],
                    "TMPFTB":[0,0.01, 10,0.6, 15,1, 25,1, 35,0],
                    "TMNFTB":[0.,0., 3.,1.]}
        self.kiosk = VariableKiosk()
        dummyday = date(2000,1,1)
        self.assim = WOFOST_Assimilation(dummyday, self.kiosk, cropdata)
        
    def runTest(self):
        from datetime import date, timedelta
        refdate = date(1999,12,31)
        # Register DVS and LAI with the kiosk
        self.kiosk.register_variable(1, "LAI", type="S", publish=True)
        self.kiosk.register_variable(1, "DVS", type="S", publish=True)
        for (IDAY,LAT,DVS,IRRAD,TMINRA,DTEMP,LAI,PGASS) in self.ref_results:
            day = refdate + timedelta(days=IDAY)
            drv = self.container()
            drv.DTEMP  = DTEMP
            drv.TMINRA = TMINRA
            drv.IRRAD  = IRRAD
            drv.LAT = LAT
            self.kiosk.set_variable(1, "LAI", LAI)
            self.kiosk.set_variable(1, "DVS", DVS)
            rpgass = self.assim(day, drv)
            self.assertAlmostEqual(rpgass, PGASS, 3)

def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_WOFOST_Assimilation))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())

