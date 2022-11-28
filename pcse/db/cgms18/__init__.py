# -*- coding: utf-8 -*-
# Copyright (c) 2004-2016 Alterra, Wageningen-UR
# Steven Hoek (steven.hoek@wur.nl), April 2016
"""Tools for reading weather data and timer, soil and site parameters
from a CGMS18 compatible database.

CGMS18 is nearly compatible with CGMS14, except for the components defined in
cgms18.dataproviders.
"""

from ..cgms14 import *
from .data_providers import SoilDataIterator
