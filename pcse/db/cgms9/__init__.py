# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Tools for reading  weather data and parameters from a CGMS9 database

Source code for the readers comes either from the PCSE data providers or the CGMS11 data providers.
"""
from ..pcse.db_input import GridWeatherDataProvider
from ..cgms11.data_providers import SoilDataProviderSingleLayer as SoilDataProvider
from ..cgms11.data_providers import SiteDataProvider, STU_Suitability, CropDataProvider
from ..cgms11.data_providers import SoilDataIterator as SoilDataIterator_CGMS11
from .data_providers import AgroManagementDataProvider


class SoilDataIterator(SoilDataIterator_CGMS11):
    """Soil data iterator for CGMS9.

    The only difference is that in CGMS9 the table is called 'ELEMENTARY_MAPPING_UNIT' and
    in CGMS11 it is called 'EMU'
    """
    emu_table_name = "elementary_mapping_unit"
