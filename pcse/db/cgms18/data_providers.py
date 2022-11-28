# -*- coding: utf-8 -*-
# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2016

"""
Data providers for weather, soil, crop, timer and site data. Also
a class for testing STU suitability for a given crop.

Data providers for CGMS18 are mostly compatible with a CGMS 14 database schema
some differences are implemented here.
"""
import datetime as dt
import os

import numpy as np
import yaml
from sqlalchemy import MetaData, select, Table, and_

from ... import exceptions as exc
from ...base import WeatherDataContainer, WeatherDataProvider
from ...util import wind10to2, safe_float, check_date, reference_ET
from ... import settings
from .. import wofost_parameters

from ..cgms14.data_providers import SoilDataIterator as SoilDataIterator14
from ..cgms14.data_providers import SoilDataProviderSingleLayer


class SoilDataIterator(SoilDataIterator14):
    """Soil data iterator for CGMS18.

    """
    tbl_link_sm_grid_cover = "link_cover_grid_smu"

