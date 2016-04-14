# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014


from sqlalchemy import create_engine, MetaData, Table

from . import db
from .models import Wofost71_WLP_FD, Wofost71_PP

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
    :param clear_table: If set to True: delete everything from the tables
        `sim_results_timeseries` and `sim_results_summary`(defaults to False)
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
    sitedata = db.pcse.fetch_sitedata(db_metadata, grid, year)
    timerdata = db.pcse.fetch_timerdata(db_metadata, grid, year, crop)
    cropdata = db.pcse.fetch_cropdata(db_metadata, grid, year, crop)
    soildata = db.pcse.fetch_soildata(db_metadata, grid)

    startdate = timerdata["START_DATE"]
    enddate = timerdata["END_DATE"]
    wdp = db.pcse.GridWeatherDataProvider(db_metadata, grid_no=grid,
                startdate=startdate, enddate=enddate)
                             
    # Initialize PCSE/WOFOST
    mode = mode.strip().lower()
    if mode == 'pp':
        wofsim = Wofost71_PP(sitedata, timerdata, soildata, cropdata, wdp)
    elif mode == 'wlp':
        wofsim = Wofost71_WLP_FD(sitedata, timerdata, soildata, cropdata, wdp)
    else:
        msg = "Unrecognized mode keyword: '%s' should be one of 'pp'|'wlp'" % mode
        raise RuntimeError(msg, mode)

    wofsim.run_till_terminate()
    
    runid = {"grid_no":grid, "crop_no":crop, "year":year, "member_id":0,
             "simulation_mode":mode}
    wofsim.store_to_database(db_metadata, runid)

