# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024

import datetime as dt

import pandas as pd
from sqlalchemy import create_engine, MetaData, Table
import sqlalchemy as sa

from ..db.pcse import fetch_cropdata, fetch_sitedata, fetch_soildata, GridWeatherDataProvider, \
    AgroManagementDataProvider
from ..base import ParameterProvider
from ..models import Wofost72_WLP_CWB, Wofost72_PP
from .. import exceptions as exc
from ..util import merge_dict


def store_to_database(engine, output, summary_output, runid):
    """Stores saved variables of the model run in a database table.

    :param meta: An SQLAlchemy metadata object providing access to the
                     database where the table 'sim_results_timeseries' can be
                     found.
    :param runid:    A dictionary providing the values for the database
                     columns that 'describe' the WOFOST run. For CGMS this
                     would be the CROP_NO, GRID_NO, YEAR thus the runid
                     would be for example be:
                     `runid={'GRID_NO':1000, 'CROP_NO':1, 'YEAR':2000}`

    Note that the records are written directly to this table. No checks on
    existing records are being carried out.
    """

    if not isinstance(runid, dict):
        msg = ("Keyword 'runid' should provide the database columns " +
               "describing the WOFOST run.")
        raise exc.PCSEError(msg)

    # Merge records with output ad summary_output variables with the run_id
    recs_output = [merge_dict(rec, runid) for rec in output]
    df_output = pd.DataFrame(recs_output)
    df_output = df_output.drop(columns=["RFTRA", "WWLOW"], errors="ignore")
    recs_summary_output = [merge_dict(rec, runid) for rec in summary_output]
    df_summary_output = pd.DataFrame(recs_summary_output)
    df_summary_output = df_summary_output.drop(columns=["CEVST", "WWLOW"], errors="ignore")

    with engine.begin() as DBconn:
        df_output.to_sql("sim_results_timeseries", DBconn, if_exists='append', index=False)
        df_summary_output.to_sql("sim_results_summary", DBconn, if_exists="append", index=False)


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
    engine = create_engine(dsn)
    meta = MetaData()
    meta.reflect(bind=engine)
    table_sim_results_ts = meta.tables['sim_results_timeseries']
    table_sim_results_smry = meta.tables['sim_results_summary']
    if clear_table is True:
        with engine.begin() as DBconn:
            DBconn.execute(table_sim_results_ts.delete())
            DBconn.execute(table_sim_results_smry.delete())
    
    # Get input data from database
    sited = fetch_sitedata(engine, meta, grid, year)
    cropd = fetch_cropdata(engine, meta, grid, year, crop)
    soild = fetch_soildata(engine, meta, grid)
    parameters = ParameterProvider(sitedata=sited, cropdata=cropd, soildata=soild)

    # Get Agromanagement
    agromanagement = AgroManagementDataProvider(engine, meta, grid, crop, year)

    # Get weather data
    wdp = GridWeatherDataProvider(engine, grid_no=grid)
                             
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
    output = wofsim.get_output()
    df = pd.DataFrame(output)
    summary_output = wofsim.get_summary_output()
    
    runid = {"grid_no":grid, "crop_no":crop, "year":year, "member_id":0,
             "simulation_mode":mode}
    store_to_database(engine, output, summary_output, runid)

