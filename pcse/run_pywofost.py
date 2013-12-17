
import sys, os
import datetime
import logging

from sqlalchemy import create_engine, MetaData, Table

from . import db
from . import util
from .pywofost import PyWofost

def run_pywofost(dsn, crop, grid, year, mode, clear_table=False,
                 multi_layer=False):
    """Provides a convenient interface for running PyWofost
    
    Starting run_pywofost() will start a PyWOFOST instance and let it
    run for 300 days for the given grid, crop, year and mode. Optionally
    it can delete everything from the table `pywofost_output`.
    
    :param dsn: PyWofost DB as SQLAlchemy data source name
    :param crop: crop number
    :param grid: grid number
    :param year: year to start
    :param mode: production mode ('pp' or 'wlp')
    :param clear_table: If set to True: delete everything from the table
        `pywofost_output` (defaults to  False)
    :param multi_layer: Use multi-layer waterbalance instead of classic
        waterbalance
    """

    # Open database connection and empty output table
    pywofost_engine = create_engine(dsn)
    pywofost_metadata = MetaData(pywofost_engine)
    table_pywofost_output = Table('pywofost_output', pywofost_metadata,
                                  autoload=True)
    if clear_table is True:
        table_pywofost_output.delete().execute()
    
    # Get input data from database
    sitedata  = db.pcse.fetch_sitedata(pywofost_metadata, grid, year)
    timerdata = db.pcse.fetch_timerdata(pywofost_metadata, grid, year, crop)
    cropdata = db.pcse.fetch_cropdata(pywofost_metadata, grid, year, crop)
    if multi_layer is False:
        soildata = db.pcse.fetch_soildata(pywofost_metadata, grid)
    else:
        soildata = db.pcse.fetch_soildata_layered(pywofost_metadata, grid)

    startdate = timerdata["START_DATE"]
    enddate = startdate + datetime.timedelta(days=365)
    meteof = db.pcse.GridWeatherDataProvider(pywofost_metadata, grid_no=grid,
                startdate=startdate, enddate=enddate)
                             
    # Run simulation
    wofsim = PyWofost(startdate, sitedata, timerdata, soildata, cropdata, meteof,
                      mode=mode, metadata=pywofost_metadata)
    
    wofsim.grow(days=365)
    
    runid = {"grid_no":grid, "crop_no":crop, "year":year, "member_id":0,
             "simulation_mode":mode}
    wofsim.store_to_database(pywofost_metadata, runid)

