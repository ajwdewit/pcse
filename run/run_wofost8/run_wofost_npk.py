# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Module defines code for running unittests on the complete PCSE/WOFOST_NPK model.
"""
import os, sys
import yaml

import pandas as pd

from pcse.models import Wofost80_NWLP_FD_beta, Wofost72_PP, Wofost80_PP_beta
from pcse.base import ParameterProvider
from pcse.fileinput import CABOFileReader, CABOWeatherDataProvider, ExcelWeatherDataProvider

test_data_dir =  os.path.join(r"C:\Users\bergh026\source\repos\pcse\run\run_wofost8", "data")

def main():
    agro = yaml.safe_load(open(os.path.join(test_data_dir, "wofost_npk.agro")))['AgroManagement']
    soil = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.soil"))
    site = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.site"))
    crop = CABOFileReader(os.path.join(test_data_dir, "wofost_npk.crop"))
    weather = CABOWeatherDataProvider("NL1", test_data_dir)

    parvalues = ParameterProvider(sitedata=site, soildata=soil, cropdata=crop)
    wofost = Wofost80_PP_beta(parvalues,  weather, agromanagement=agro)
    wofost.run_till_terminate()
    output_potential = pd.DataFrame(wofost.get_output()).set_index("day")
    output_potential.to_excel(os.path.join("c:/temp", "WOFOST80_PP.xlsx"))

    wofost = Wofost80_NWLP_FD_beta(parvalues,  weather, agromanagement=agro)
    wofost.run_till_terminate()
    output_limited = pd.DataFrame(wofost.get_output()).set_index("day")
    output_limited.to_excel(os.path.join("c:/temp", "WOFOST80_NWLP.xlsx"))

    wofost = Wofost72_PP(parvalues,  weather, agromanagement=agro)
    wofost.run_till_terminate()
    output_potential = pd.DataFrame(wofost.get_output()).set_index("day")
    output_potential.to_excel(os.path.join("c:/temp", "WOFOST72_PP.xlsx"))

if __name__ == '__main__':
   main()
