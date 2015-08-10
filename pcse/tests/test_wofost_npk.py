# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Module defines code for running unittests on the complete PCSE/WOFOST_NPK model.
"""
import os, sys
import unittest
import csv
import yaml

from ..engine import Engine
from ..base_classes import MultiCropParameterProvider
from ..fileinput import CABOFileReader, CABOWeatherDataProvider

test_data_dir =  os.path.join(os.path.dirname(__file__), "test_data")


class TestWOFOSTNPK_WinterWheat(unittest.TestCase):

    def setUp(self):
        amgt = yaml.load(open(os.path.join(test_data_dir, "wofost_npk.amgt")))['AgroManagement']
        soil = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.soil"))
        site = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.site"))
        crop = {"winter-wheat": CABOFileReader(os.path.join(test_data_dir, "wofost_npk.crop"))}
        weather = CABOWeatherDataProvider("NL1", test_data_dir)

        parvalues = MultiCropParameterProvider(sitedata=site, soildata=soil, multi_cropdata=crop)
        wofost = Engine(parvalues,  weather, agromanagement=amgt, config="Wofost71_NPK.conf")
        wofost.run(days=300)
        self.output = wofost.get_output()

    def _safe_cast_float(self, ref):
        """The reference data come in as strings from the CSV file.
        This function tries to convert to float or returns None
        """
        try:
            return float(ref)
        except ValueError as exc:
            return None

    def runTest(self):
        refdata_file = os.path.join(test_data_dir, "wofost_npk_reference_results.csv")
        ref_results = []
        with open(refdata_file, "rb") as fp:
            reader = csv.DictReader(fp)
            ref_results = [row for row in reader]
        msg = "Different number of rows in model output and reference results"
        self.assertEqual(len(ref_results), len(self.output), msg)
        ntests = 0
        for ref, out in zip(ref_results, self.output):
            msg = "Date not equal when comparing reference with simulation."
            day = ref.pop("day")
            self.assertEqual(day, str(out.pop("day")), msg)
            for name, refvalue in ref.items():
                rv = self._safe_cast_float(refvalue)
                sv = out[name]
                if rv is None and sv is None:
                    continue
                msg = "%s, %s: %f != %f" % (day, name, sv, rv)
                self.assertAlmostEqual(rv, sv, msg=msg, delta=0.001)
            ntests += 1
        msg = "Number of days tested should be 216, found %i" % ntests
        self.assertEqual(ntests, 216, msg)



def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestWOFOSTNPK_WinterWheat))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())
