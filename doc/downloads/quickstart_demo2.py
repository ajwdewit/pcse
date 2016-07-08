# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), June 2017

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pcse
from pcse.db import NASAPowerWeatherDataProvider
from pcse.fileinput import CABOFileReader
from pcse.fileinput import PCSEFileReader
from pcse.fileinput import YAMLAgroManagementReader
from pcse.models import Wofost71_WLP_FD
from pcse.base_classes import ParameterProvider

# First set the location where the crop, soil and crop calendar files can be found
data_dir = r""
if data_dir == "":
    print("Variable 'data_dir' in line 19 must be set to the location of the the data folder")
    sys.exit()

# Retrieve weather data from the NASA Power database
wdp = NASAPowerWeatherDataProvider(latitude=52, longitude=5)

# Read parameter values from the input files
cropdata = CABOFileReader(os.path.join(data_dir,'sug0601.crop'))
soildata = CABOFileReader(os.path.join(data_dir,'ec3.soil'))
sitedata = {'SSMAX'  : 0.,
            'IFUNRN' : 0,
            'NOTINF' : 0,
            'SSI'    : 0,
            'WAV'    : 100,
            'SMLIM'  : 0.03,
            'CO2'    : 360}
parameters = ParameterProvider(cropdata=cropdata, soildata=soildata, sitedata=sitedata)

# Read agromanagement
agromanagement = YAMLAgroManagementReader(os.path.join(data_dir,'sugarbeet_calendar.amgt'))

# Start WOFOST
wf = Wofost71_WLP_FD(parameters, wdp, agromanagement)
wf.run_till_terminate()

# Get time-series output from WOFOST and take the selected variables
output = wf.get_output()
varnames = ["day", "DVS", "TAGP", "LAI", "SM"]
tmp = {}
for var in varnames:
    tmp[var] = [t[var] for t in output]
day = tmp.pop("day")

# make a figure with 2x2 subplots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,8))
# loop over variables and axes in order to plot on each axes
for var, ax in zip(["DVS", "TAGP", "LAI", "SM"], axes.flatten()):
    ax.plot_date(day, tmp[var], 'b-')
    ax.set_title(var)
# Autoformat the dates on the x-axis and generate a .PNG file
fig.autofmt_xdate()
fig.savefig('sugarbeet.png')

