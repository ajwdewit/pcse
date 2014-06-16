# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), May 2014
"""Tools for reading  weather files in the NetCDF4 format used
for crop simulation models:
- NetcdfWeatherDataProvider reads NETCDF4 files for use in PyWOFOST
"""
from .netcdf4_reader import NetcdfWeatherDataProvider;