# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

import datetime as dt
from sqlalchemy import create_engine, MetaData, Table
import sqlalchemy as sa

from ..db.pcse import fetch_cropdata, fetch_sitedata, fetch_soildata, GridWeatherDataProvider, \
    AgroManagementDataProvider
from ..base_classes import ParameterProvider
from ..models import Wofost71_WLP_FD, Wofost71_PP
from .. import exceptions as exc
from ..util import merge_dict


def store_to_database(metadata, output, summary_output, runid):
    """Stores saved variables of the model run in a database table.

    :param metadata: An SQLAlchemy metadata object providing access to the
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

    if not isinstance(metadata, sa.schema.MetaData):
        msg = ("Keyword metadata should provide an SQLAlchemy " +
               "MetaData object.")
        raise exc.PCSEError(msg)

    # Merge records with output ad summary_output variables with the run_id
    recs_output = [merge_dict(rec, runid) for rec in output]
    recs_summary_output = [merge_dict(rec, runid) for rec in summary_output]

    table_sim_results_ts = sa.Table('sim_results_timeseries', metadata,
                                    autoload=True)
    i = table_sim_results_ts.insert()
    i.execute(recs_output)

    table_sim_results_smry = sa.Table('sim_results_summary', metadata,
                                      autoload=True)
    i = table_sim_results_smry.insert()
    i.execute(recs_summary_output)


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
    db_engine = create_engine(dsn)
    db_metadata = MetaData(db_engine)
    table_sim_results_ts = Table('sim_results_timeseries', db_metadata,
                                 autoload=True)
    table_sim_results_smry = Table('sim_results_summary', db_metadata,
                                   autoload=True)
    if clear_table is True:
        table_sim_results_ts.delete().execute()
        table_sim_results_smry.delete().execute()
    
    # Get input data from database
    sited = fetch_sitedata(db_metadata, grid, year)
    cropd = fetch_cropdata(db_metadata, grid, year, crop)
    soild = fetch_soildata(db_metadata, grid)
    parameters = ParameterProvider(sitedata=sited, cropdata=cropd, soildata=soild)

    # Get Agromanagement
    agromanagement = AgroManagementDataProvider(db_engine, grid, crop, year)

    start_date = agromanagement[0].keys()[0]
    end_date = start_date + dt.timedelta(days=365)
    wdp = GridWeatherDataProvider(db_engine, grid_no=grid, start_date=start_date, end_date=end_date)
                             
    # Initialize PCSE/WOFOST
    mode = mode.strip().lower()
    if mode == 'pp':
        wofsim = Wofost71_PP(parameters, wdp, agromanagement)
    elif mode == 'wlp':
        wofsim = Wofost71_WLP_FD(parameters, wdp, agromanagement)
    else:
        msg = "Unrecognized mode keyword: '%s' should be one of 'pp'|'wlp'" % mode
        raise RuntimeError(msg, mode)

    wofsim.run_till_terminate()
    output = wofsim.get_output()
    summary_output = wofsim.get_summary_output()
    
    runid = {"grid_no":grid, "crop_no":crop, "year":year, "member_id":0,
             "simulation_mode":mode}
    store_to_database(db_metadata, output, summary_output, runid)

