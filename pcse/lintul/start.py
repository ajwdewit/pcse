# -*- coding: utf-8 -*-
from pcse.base_classes import ParameterProvider
from pcse.engine import Engine
from pcse.fileinput.cabo_weather import CABOWeatherDataProvider
from datetime import date
import lintul3parameters
from pcse.lintul.lintul3 import SubModel
from numbers import Number


class Lintul3Model(Engine):

    

    @classmethod
    def readModelParameters(cls, module):
        '''
        read model parameters from a slightly (syntax) adapted Model.dat
        :param module: an imported module
        '''
        allvars = module.__dict__
        params = dict((k, v) for k, v in allvars.items() if  not k.startswith('_'))
        return params        


    @classmethod
    def start(cls, year=1997, outputProc = None):
        
        # Timer: starting day, final day and model output
        timerdata           = {
                               "CAMPAIGNYEAR": year,                      # year of the agricultural campaign (e.g. harvest year)
                                 "START_DATE": date(year, 01, 01),        # date of the start of the simulation
                                   "END_DATE": date(year, 12, 31),        # date last possible day of the simulation
                            "CROP_START_TYPE": "emergence",               # 'emergence' or 'sowing'
                            "CROP_START_DATE": date(year, 03, 31),        # date of the start of the crop simulation
                              "CROP_END_TYPE": "earliest",                # 'maturity' | 'harvest' |'earliest'
                              "CROP_END_DATE": date(year, 10, 20),        # date of the end of the crop simulation in case of CROP_END_TYPE == 'harvest' | 'earliest'
                               "MAX_DURATION": 366                        # maximum number of days of the crop simulation
                               }
        parameterprovider   = ParameterProvider(Lintul3Model.readModelParameters(lintul3parameters), timerdata, {}, {})
        weatherdataprovider = CABOWeatherDataProvider("NL1", "D:/Projects/pcse/lintul/Lintul-3 model/data/", ETmodel='P')
    
        SubModel.onOutput   = outputProc    
        return Lintul3Model(parameterprovider, weatherdataprovider, config="lintul3.conf.py")

    

if (__name__ == "__main__"):
    
    class P:            
        __lineBuffer = {}
        __headerBuffer = {}
        __headerPrinted = False


        def printRow(self, values):
            for v in values:
                if isinstance(v, Number):
                    print "%f\t" % (v),
                else:
                    print "%s\t" % (v),
            print
            
            
        def __call__(self, values, header = None):
            if header:
                self.printRow(header)
            self.printRow(values)


        
    p = P()
    
    sim = Lintul3Model.start(1987, outputProc=p)
    
    sim.run(365)
    SubModel.doOutput(sim.crop, 99999, {})
    print "END STOP ENDJOB "*5

    
