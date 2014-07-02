# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), May 2014
"""Tools for reading  weather files in the NetCDF4 format used
for crop simulation models:
- NetcdfWeatherDataConverter converts NETCDF4 files to HDF5 files for use in PyWOFOST
"""
from .netcdf4converter import NetcdfWeatherDataConverter;