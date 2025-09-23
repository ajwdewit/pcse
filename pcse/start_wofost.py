# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
import sys, os
from collections import namedtuple
import sqlite3

from .tests.db_input import GridWeatherDataProvider, fetch_soildata, fetch_sitedata, fetch_cropdata, \
    AgroManagementDataProvider
from .base import ParameterProvider
from .models import Wofost72_PP, Wofost72_WLP_CWB
from .settings import settings


def namedtuple_factory(cursor, row):
    """Creates an SQLite row factor for named tuples.

    see: https://docs.python.org/3/library/sqlite3.html#how-to-create-and-use-row-factories
    """
    fields = [column[0] for column in cursor.description]
    cls = namedtuple("Row", fields)
    return cls._make(row)


def start_wofost(grid=31031, crop=1, year=2000, mode='wlp'):
    """Provides a convenient interface for starting a WOFOST instance for the internal Demo DB.
    
    If started with no arguments, the routine will connnect to the
    demo database and initialize WOFOST for winter-wheat (cropno=1)
    in Spain (grid_no=31031) for the year 2000 in water-limited
    production (mode='wlp')
    
    
    :param grid: grid number, defaults to 31031
    :param crop: crop number, defaults to 1 (winter-wheat in the demo database)
    :param year: year to start, defaults to 2000
    :param mode: production mode ('pp' or 'wlp'), defaults to 'wlp'

    example::
    
        >>> import pcse
        >>> wofsim = pcse.start_wofost(grid=31031, crop=1, year=2000, 
        ...   mode='wlp')
        >>> 
        >>> wofsim
        <pcse.models.Wofost71_WLP_FD at 0x35f2550>
        >>> wofsim.run(days=300)
        >>> wofsim.get_variable('tagp')
        15261.752187075261
    """

    # Open database connections
    db_location = os.path.join(settings.PCSE_USER_HOME, "pcse.db")
    DBconn = sqlite3.connect(db_location)
    DBconn.row_factory = namedtuple_factory

    # Get input data from database
    agromanagement = AgroManagementDataProvider(DBconn, grid, crop, year)
    sited  = fetch_sitedata(DBconn, grid, year)
    cropd = fetch_cropdata(DBconn, grid, year, crop)
    soild = fetch_soildata(DBconn, grid)
    parvalues = ParameterProvider(sitedata=sited, soildata=soild, cropdata=cropd)

    wdp = GridWeatherDataProvider(DBconn, grid_no=grid)
                             
    # Initialize PCSE/WOFOST
    mode = mode.strip().lower()
    if mode == 'pp':
        wofsim = Wofost72_PP(parvalues, wdp, agromanagement)
    elif mode == 'wlp':
        wofsim = Wofost72_WLP_CWB(parvalues, wdp, agromanagement)
    else:
        msg = "Unrecognized mode keyword: '%s' should be one of 'pp'|'wlp'" % mode
        raise RuntimeError(msg)
    return wofsim
