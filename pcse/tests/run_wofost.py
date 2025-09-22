# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024

import sqlite3
from collections import namedtuple
import pandas as pd

from .db_input import fetch_cropdata, fetch_sitedata, fetch_soildata, AgroManagementDataProvider, GridWeatherDataProvider
from ..base import ParameterProvider
from ..models import Wofost72_WLP_CWB, Wofost72_PP


def namedtuple_factory(cursor, row):
    """Creates an SQLite row factor for named tuples.

    see: https://docs.python.org/3/library/sqlite3.html#how-to-create-and-use-row-factories
    """
    fields = [column[0] for column in cursor.description]
    cls = namedtuple("Row", fields)
    return cls._make(row)


def run_wofost(dsn, crop, grid, year, mode, clear_table=False):
    """Provides a convenient interface for running PCSE/WOFOST from a
    PCSE database.
    
    Starting run_wofost() will start a PCSE/WOFOST instance and let it
    run until it terminates for the given grid, crop, year and mode.
    Optionally it can delete everything from tables `sim_results_timeseries`
    and `sim_results_summary`.
    
    :param dsn: PCSE DB as SQLAlchemy data source name
    :param crop: crop number
    :param grid: grid number
    :param year: year to start
    :param mode: production mode ('pp' or 'wlp')
    :param clear_table: If set to True delete everything from the tables
        'sim_results_timeseries' and  'sim_results_summary' (defaults to False)
    """

    # Open database connection and empty output table
    DBconn = sqlite3.connect(dsn)
    DBconn.row_factory = namedtuple_factory
    cursor = DBconn.cursor()
    if clear_table is True:
        cursor.execute("delete from sim_results_timeseries")
        cursor.execute("delete from sim_results_summary")
        cursor.close()

    # Get input data from database
    sited = fetch_sitedata(DBconn, grid, year)
    cropd = fetch_cropdata(DBconn, grid, year, crop)
    soild = fetch_soildata(DBconn, grid)
    parameters = ParameterProvider(sitedata=sited, cropdata=cropd, soildata=soild)

    # Get Agromanagement
    agromanagement = AgroManagementDataProvider(DBconn, grid, crop, year)

    # Get weather data
    wdp = GridWeatherDataProvider(DBconn, grid_no=grid)
                             
    # Initialize PCSE/WOFOST
    mode = mode.strip().lower()
    if mode == 'pp':
        wofsim = Wofost72_PP(parameters, wdp, agromanagement)
    elif mode == 'wlp':
        wofsim = Wofost72_WLP_CWB(parameters, wdp, agromanagement)
    else:
        msg = "Unrecognized mode keyword: '%s' should be one of 'pp'|'wlp'" % mode
        raise RuntimeError(msg, mode)

    wofsim.run_till_terminate()
    df_output = pd.DataFrame(wofsim.get_output())
    df_summary_output = pd.DataFrame(wofsim.get_summary_output())
    
    runid = {"grid_no":grid, "crop_no":crop, "year":year, "member_id":0,
             "simulation_mode":mode}
    for name, value in runid.items():
        df_output[name] = value
        df_summary_output[name] = value

    columns = ["grid_no", "crop_no", "year", "day", "simulation_mode", "member_id",
               "DVS", "LAI", "TAGP", "TWSO", "TWLV", "TWST", "TWRT", "TRA", "RD", "SM"]
    df_output = df_output[columns]
    df_output.to_sql("sim_results_timeseries", DBconn, index=False, if_exists='append')

    columns = ["grid_no", "crop_no", "year", "simulation_mode", "member_id",
               "DVS", "LAIMAX", "TAGP", "TWSO", "TWLV", "TWST", "TWRT", "CTRAT",
               "RD", "DOS", "DOE", "DOA", "DOM", "DOH", "DOV"]
    df_summary_output = df_summary_output[columns]
    df_summary_output.to_sql("sim_results_summary", DBconn, index=False, if_exists='append')

