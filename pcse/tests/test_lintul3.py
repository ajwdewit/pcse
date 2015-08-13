'''
Created on 30 Apr 2015

@author: winte005
'''
import unittest

from pcse.lintul.start import Lintul3Model
import os.path


class Test(unittest.TestCase):


    test_input = "lintul3.csv"
    
    def setUp(self):
        test_dir = os.path.dirname(__file__)
        csv_file = os.path.join(test_dir, "test_data", self.test_input)
        self.csv_data = {}
        header = None
        for s in open(csv_file):
            fields = s.split(',')
            if (header == None):
                header = fields
            else:
                # convert numerals:
                for i in range(len(fields)):
                    try:
                        fields[i] = int(fields[i])
                    except:
                        try:
                            fields[i] = float(fields[i])
                        except:
                            pass
                        
                time = fields[0]
                self.csv_data[time] = dict(zip(header, fields))
                


    class CompareOutput:            
        __lineBuffer = {}
        __headerBuffer = {}
        __headerPrinted = False

        def __init__(self, owner):
            self.owner = owner

            
            
        def __call__(self, values, header = None):
            if header:
                self.header = header

            def sameValue(a, b, rtol = 0., atol = 0.):
                # contrary to Numpy.allequal this function is symmetric for a and b
                equality = abs(abs(a) - abs(b)) <= (atol + rtol * max(abs(b), abs(a)))
                return equality
                # -------------
            n = 0
            hdr = self.header
            assert len(values) == len(hdr), "row has %d elements and header has %d" %(len(values), len(hdr))
            time = values[0]
            data = self.owner.csv_data[time]
            for i in range(len(values)):
                key = hdr[i]
                print key,
                if data.has_key(key):
                    # xBALAN balances suffer from rounding errors... 
                    rtol = 0.001 
                    atol = 0.001 if key.endswith("BALAN") else 0.000
                    equality = sameValue(values[i], data[key], rtol, atol)
                                
                    assert  equality, "values differ at time = %d for %s... found: %f expected: %f"  % (time, key, values[i], data[key])
                    n += 1
                    
            print "%d values checked at time %d" %(n, time)
        
        
    # obsolete
    def compareOutput(self, aRow, header = None):

        def sameValue(a, b, rtol = 0., atol = 0.):
            # contrary to Numpy.allequal this function is symmetric for a and b
            equality = abs(abs(a) - abs(b)) <= (atol + rtol * max(abs(b), abs(a)))
            return equality
            # -------------
        
        if header:
            self.header = header

        hdr = self.header
        assert len(aRow) == len(hdr), "row has %d elements and header has %d" %(len(aRow), len(hdr))
        time = aRow[0]
        data = self.csv_data[time]
        n = 0
        for i in range(len(aRow)):
            key = hdr[i]
            if data.has_key(key) and not (aRow[i] == '-'):
                # xBALAN balances suffer from rounding errors... 
                rtol = 0.001 
                atol = 0.001 if key.endswith("BALAN") else 0.000
                equality = sameValue(aRow[i], data[key], rtol, atol)
                            
                assert  equality, "values differ at time = %d for %s... found: %f expected: %f"  % (time, key, aRow[i], data[key])
                n+=1
                
        print "%d values checked at time %d" %(n, time)
        
        
        
    def testName(self):
        callBack = self.CompareOutput(self)
        sim = Lintul3Model.start(year=1987, outputProc=callBack)
        
        sim.run(182)



    
    
if __name__ == "__main__":
    unittest.main()