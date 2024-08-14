# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
"""Tools for reading  weather and parameter files.

Contains the following data providers:
- CABOWeatherDataProvider reads CABOWE weather files for use in PCSE
- CABOFileReader reads CABO parameter files.
- PCSEFileReader reads parameters files in the PCSE format
- ExcelWeatherDataProvider reads weather data in xlsx format
- CSVWeatherDataProvider reads weather data in CSV format
- YAMLAgroManagementReader reads agromanagement data in YAML format
- YAMLCropDataProvider reads crop parameters in YAML format

Site data providers for several WOFOST versions:
- WOFOST72SiteDataProvider
- WOFOST73SiteDataProvider
- WOFOST81SiteDataProvider_Classic
- WOFOST81SiteDataProvider_SNOMIN


"""

from .cabo_reader import CABOFileReader
from .cabo_weather import CABOWeatherDataProvider
from .pcsefilereader import PCSEFileReader
from .excelweatherdataprovider import ExcelWeatherDataProvider
from .csvweatherdataprovider import CSVWeatherDataProvider
from .yaml_agro_loader import YAMLAgroManagementReader
from .yaml_cropdataprovider import YAMLCropDataProvider
from .nasapower import NASAPowerWeatherDataProvider
from .sitedataproviders import WOFOST72SiteDataProvider, WOFOST73SiteDataProvider, \
    WOFOST81SiteDataProvider_Classic, WOFOST81SiteDataProvider_SNOMIN
from .soildataproviders import DummySoilDataProvider