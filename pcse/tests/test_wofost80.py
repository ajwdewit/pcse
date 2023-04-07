# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Module defines code for running unittests on the complete PCSE/WOFOST_NPK model.
"""
import os, sys
import unittest
import csv
import yaml

from ..models import Wofost80_PP_beta, Wofost80_NWLP_FD_beta
from ..base import ParameterProvider
from ..fileinput import CABOFileReader, CABOWeatherDataProvider

test_data_dir =  os.path.join(os.path.dirname(__file__), "test_data")


class TestWOFOST80_Potential_WinterWheat(unittest.TestCase):
    reference_results = "wofost_npk_potential.csv"
    model = Wofost80_PP_beta

    def setUp(self):
        agro = yaml.safe_load(open(os.path.join(test_data_dir, "wofost_npk.agro")))['AgroManagement']
        soil = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.soil"))
        site = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.site"))
        crop = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.crop"))
        weather = CABOWeatherDataProvider("NL1", test_data_dir)

        parvalues = ParameterProvider(sitedata=site, soildata=soil, cropdata=crop)
        wofost = self.model(parvalues,  weather, agromanagement=agro)
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
        refdata_file = os.path.join(test_data_dir, self.reference_results)
        with open(refdata_file, "r") as fp:
            reader = csv.DictReader(fp)
            ref_results = [row for row in reader]
        msg = "Different number of rows in model output and reference results"
        self.assertEqual(len(ref_results), len(self.output), msg)
        ntests = 0
        for ref_results, sim_outputs in zip(ref_results, self.output):
            msg = "Date not equal when comparing reference with simulation."
            day = ref_results.pop("day")
            self.assertEqual(day, str(sim_outputs.pop("day")), msg)
            for name, ref_value in ref_results.items():
                if name not in sim_outputs:
                    continue
                v1 = self._safe_cast_float(ref_value)
                v2 = sim_outputs[name]
                if v1 is None and v2 is None:
                    continue
                msg = "%s, %s: %f != %f" % (day, name, v2, v1)
                self.assertAlmostEqual(v1, v2, msg=msg, delta=0.001)
            ntests += 1
        msg = "Number of days tested should be 216, found %i" % ntests
        self.assertEqual(ntests, 216, msg)


class TestWOFOST80_WaterLimited_WinterWheat(TestWOFOST80_Potential_WinterWheat):
    reference_results = "wofost_npk_waterlimited.csv"
    model = Wofost80_NWLP_FD_beta


def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestWOFOST80_Potential_WinterWheat))
    suite.addTest(unittest.makeSuite(TestWOFOST80_WaterLimited_WinterWheat))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())
