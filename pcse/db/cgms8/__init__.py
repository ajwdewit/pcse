# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Tools for reading  weather data and parameters from a CGMS8 database

Source code partially comes from the CGMS12 data providers as the database has a similar structure
for many of the soil and crop related tables.
"""
from ..cgms12.data_providers import SoilDataProviderSingleLayer as SoilDataProvider
from ..cgms12.data_providers import SiteDataProvider, STU_Suitability, CropDataProvider
from ..cgms12.data_providers import SoilDataIterator as SoilDataIterator_CGMS12
from .data_providers import AgroManagementDataProvider, GridWeatherDataProvider


class SoilDataIterator(SoilDataIterator_CGMS12):
    """Soil data iterator for CGMS8.

    The only difference is that in CGMS8 the table is called 'ELEMENTARY_MAPPING_UNIT' and
    in CGMS12 it is called 'EMU'
    """
    emu_table_name = "elementary_mapping_unit"
