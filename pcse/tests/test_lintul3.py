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

test_data_dir =  os.path.join(os.path.dirname(__file__), "test_data")

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
        self.assertEqual(ntests, 225, msg)



def suite():
    """ This defines all the tests of a module"""
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestLINTUL3_SpringWheat))
    return suite

if __name__ == '__main__':
   unittest.TextTestRunner(verbosity=2).run(suite())

#
#
# class Test(unittest.TestCase):
#
#
#     test_input = "lintul3.csv"
#
#     def setUp(self):
#         test_dir = os.path.dirname(__file__)
#         csv_file = os.path.join(test_dir, "test_data", self.test_input)
#         self.csv_data = {}
#         header = None
#         for s in open(csv_file):
#             fields = s.split(',')
#             if (header == None):
#                 header = fields
#             else:
#                 # convert numerals:
#                 for i in range(len(fields)):
#                     try:
#                         fields[i] = int(fields[i])
#                     except:
#                         try:
#                             fields[i] = float(fields[i])
#                         except:
#                             pass
#
#                 time = fields[0]
#                 self.csv_data[time] = dict(zip(header, fields))
#
#
#
#     class CompareOutput:
#         __lineBuffer = {}
#         __headerBuffer = {}
#         __headerPrinted = False
#
#         def __init__(self, owner):
#             self.owner = owner
#
#
#
#         def __call__(self, values, header = None):
#             if header:
#                 self.header = header
#
#             def sameValue(a, b, rtol = 0., atol = 0.):
#                 # contrary to Numpy.allequal this function is symmetric for a and b
#                 equality = abs(abs(a) - abs(b)) <= (atol + rtol * max(abs(b), abs(a)))
#                 return equality
#                 # -------------
#             n = 0
#             hdr = self.header
#             assert len(values) == len(hdr), "row has %d elements and header has %d" %(len(values), len(hdr))
#             time = values[0]
#             data = self.owner.csv_data[time]
#             for i in range(len(values)):
#                 key = hdr[i]
#                 print key,
#                 if data.has_key(key):
#                     # xBALAN balances suffer from rounding errors...
#                     rtol = 0.001
#                     atol = 0.001 if key.endswith("BALAN") else 0.000
#                     equality = sameValue(values[i], data[key], rtol, atol)
#
#                     assert  equality, "values differ at time = %d for %s... found: %f expected: %f"  % (time, key, values[i], data[key])
#                     n += 1
#
#             print "%d values checked at time %d" %(n, time)
#
#
#     # obsolete
#     def compareOutput(self, aRow, header = None):
#
#         def sameValue(a, b, rtol = 0., atol = 0.):
#             # contrary to Numpy.allequal this function is symmetric for a and b
#             equality = abs(abs(a) - abs(b)) <= (atol + rtol * max(abs(b), abs(a)))
#             return equality
#             # -------------
#
#         if header:
#             self.header = header
#
#         hdr = self.header
#         assert len(aRow) == len(hdr), "row has %d elements and header has %d" %(len(aRow), len(hdr))
#         time = aRow[0]
#         data = self.csv_data[time]
#         n = 0
#         for i in range(len(aRow)):
#             key = hdr[i]
#             if data.has_key(key) and not (aRow[i] == '-'):
#                 # xBALAN balances suffer from rounding errors...
#                 rtol = 0.001
#                 atol = 0.001 if key.endswith("BALAN") else 0.000
#                 equality = sameValue(aRow[i], data[key], rtol, atol)
#
#                 assert  equality, "values differ at time = %d for %s... found: %f expected: %f"  % (time, key, aRow[i], data[key])
#                 n+=1
#
#         print "%d values checked at time %d" %(n, time)
#
#
#
#     def testName(self):
#         callBack = self.CompareOutput(self)
#         sim = Lintul3Model.start(year=1987, outputProc=callBack)
#
#         sim.run(182)



    
    
if __name__ == "__main__":
    unittest.main()