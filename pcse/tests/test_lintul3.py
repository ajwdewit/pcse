'''
Created on 30 Apr 2015

@author: winte005
'''
from __future__ import print_function
import os, sys
import unittest
import csv
import yaml
import pdb

from ..engine import Engine
from ..base_classes import MultiCropParameterProvider
from ..fileinput import PCSEFileReader, CABOWeatherDataProvider

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")

class TestLINTUL3_SpringWheat(unittest.TestCase):

    def setUp(self):
        amgt = yaml.load(open(os.path.join(test_data_dir, "lintul3_parameters.amgt")))['AgroManagement']
        soil = PCSEFileReader(os.path.join(test_data_dir, "lintul3_parameters.soil"))
        site = {}
        crop = {"spring-wheat": PCSEFileReader(os.path.join(test_data_dir, "lintul3_parameters.crop"))}
        weather = CABOWeatherDataProvider("NL1", test_data_dir)

        parvalues = MultiCropParameterProvider(sitedata=site, soildata=soil, multi_cropdata=crop,
                                               init_root_depth_name="ROOTDI", max_root_depth_name="ROOTDM")
        lintul3 = Engine(parvalues,  weather, agromanagement=amgt, config="lintul3.conf.py")
        lintul3.run(days=300)
        self.output = lintul3.get_output()

    def _safe_cast_float(self, ref):
        """The reference data come in as strings from the CSV file.
        This function tries to convert to float or returns None
        """
        try:
            return float(ref)
        except ValueError as exc:
            return None

    def runTest(self):
        refdata_file = os.path.join(test_data_dir, "LINTUL3_reference_results.csv")
        ref_results = []
        with open(refdata_file, "rb") as fp:
            reader = csv.DictReader(fp)
            ref_results = [row for row in reader]
        msg = "Different number of rows in model output and reference results"
        self.assertEqual(len(ref_results), len(self.output), msg)
        ntests = 0
        for ref, out in zip(ref_results, self.output):
            day = ref.pop("day")
            msg = "Date (%s) not equal when comparing reference with simulation." % day
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
        self.assertEqual(ntests, 225, msg)


def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestLINTUL3_SpringWheat))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())

