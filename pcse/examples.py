# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Examples for running PCSE/WOFOST in batch mode

Function defined here:
* run_example_simple
* run_example_ensemble
* run_example_enkf
"""
import sys, os
import datetime
import time
import socket
import logging
import logging.config

from sqlalchemy import create_engine, MetaData, Table
from numpy import random

from .run_wofost import run_wofost

def run_example1(dsn=None):
    """Simple example for running PyWOFOST
    
    Runs PyWOFOST for 6 crop types for one location in water-limited and
    potential production mode. Output is writted to the demo database and to
    a file with a name <grid>_<crop>_<year>_<mode>.out
    
    Parameters:
    dsn - SQLAlchemy data source name pointing to the database to be used.
    """

    installdir = os.path.dirname(os.path.abspath(__file__))

    logconfig = os.path.join(installdir, 'logging.conf') 
    logging.config.fileConfig(logconfig)

    if dsn is None: # Assume SQlite demo DB
        db_location = os.path.join(installdir, "database","pywofost.db")
        dsn = "sqlite:///" + db_location
    
    # Open database connection and empty output table
    pywofost_engine = create_engine(dsn)
    pywofost_metadata = MetaData(pywofost_engine)
    table_pywofost_output = Table('pywofost_output', pywofost_metadata,
                                  autoload=True)
    table_pywofost_output.delete().execute()
    
    # Loop over years, grids and crops
    years = [2000,]
    crops = [1,2,3,7,10,11]
    grids = [31031,]
    t1 = time.time()
    for year in years:
      for grid in grids:
         for crop in crops:
           print "##############################################"
           print "SIMULATING GRID: %s, CROP: %s, YEAR: %s" % (grid, crop, year)
           
           # Simulate potential production
           run_wofost(dsn, crop, grid, year, mode='PP')
           
           # Simulate water-limited production
           run_wofost(dsn, crop, grid, year, mode='WLP')
           
    logging.shutdown()
    print "Time elapsed: %7.1f seconds" % (time.time() - t1)


#-------------------------------------------------------------------------------
def initial_soil_moisture(ensemble_size, rv, SMW, SMFCF, SM0, RDMSOL):
    "Generates initial values for initial root zone soil moisture (WAV)"
    
    SM_range = (SMFCF - SMW) * RDMSOL
    init_WAV = rv.normal(SM_range/2., SM_range*0.2, (ensemble_size,))
    init_WAV = init_WAV.clip(min=0, max=SM_range)
    Esitedata = [{"WAV": WAV} for WAV in init_WAV]
    return Esitedata
#-------------------------------------------------------------------------------
def run_example_ensemble(dsn=None):
    """Example for running PyWOFOST in ensemble mode
    
    Runs PyWOFOST for 6 crop types for one location in water-limited mode
    using an ensemble of 50 members. Each member is initialized with different
    values for the initial soil moisture and each member receives different
    values for rainfall as forcing variable. Executing PyWOFOST runs is done
    through the task manager. Output is writted to the database 
    
    Parameters:
    dsn - SQLAlchemy data source name pointing to the database to be used.
    """

    from .pywofost_ensemble import PyWofostEnsemble
    from . import db_input as dbi
    from .taskmanager import TaskManager
    from sqlalchemy.exceptions import SQLAlchemyError
    
    installdir = os.path.dirname(os.path.abspath(__file__))
    logconfig = os.path.join(installdir, 'logging.conf') 
    logging.config.fileConfig(logconfig)

    if dsn is None: # Assume SQlite demo DB
        db_location = os.path.join(installdir, "database","pywofost.db")
        dsn = "sqlite:///" + db_location
    
    # Open database connection and empty output table
    engine = create_engine(dsn)
    connection = engine.connect()
    metadata = MetaData(engine)
    table_pywofost_output = Table('pywofost_output', metadata,
                                  autoload=True)
    table_pywofost_output.delete().execute()
    table_tasklist = Table("tasklist", metadata, autoload=True)
    table_tasklist.update().execute(status='Pending', hostname='')
    
    # Define ensemble size
    ensemble_size = 50

    # Initialise task manager
    taskmanager = TaskManager(metadata, connection, dbtype="SQLite",
                              hostname=socket.gethostname())
    # Loop until no tasks are left
    task = taskmanager.get_task()
    while task is not None:
        try:
            print "Running task: %i" % (task["task_id"])
            task_id = task["task_id"]
            grid_no = task["grid_no"]
            year = task["year"]
            crop_no = task["crop_no"]
            mode = task["sim_mode"]

            # Define random variable, take seed from database
            rv = random.RandomState(seed=int(task["randomseed"]))
            
            # Get Timer settings, Site, soil and Crop parameters
            timerdata = dbi.fetch_timerdata(metadata, grid_no, year,
                                                 crop_no)
            sitedata  = dbi.fetch_sitedata(metadata, grid_no, year)
            cropdata = dbi.fetch_cropdata(metadata, grid_no, year,
                                               crop_no)
            soildata = dbi.fetch_soiltype(metadata, grid_no)
            
            # Get Meteo input data
            meteo_start_date = timerdata["WB_START_DATE"]
            meteo_end_date = timerdata["CROP_END_DATE"]
            meteof = dbi.MeteoFetcher(metadata, grid_no=grid_no,
                     startdate=meteo_start_date, enddate=meteo_end_date,
                     ensemble_mode=True)
            
            # Generate ensemble of initial soil moisture values
            Esitedata = initial_soil_moisture(ensemble_size, rv, 
                        soildata["SMW"], soildata["SMFCF"],
                        soildata["SM0"], soildata["RDMSOL"])

            # Initialize ensemble
            ensemble = PyWofostEnsemble(sitedata, timerdata, soildata, cropdata,
                                        meteof, ensemble_size=ensemble_size,
                                        mode=mode, Esitedata=Esitedata,
                                        metadata=metadata)

            # Run ensemble simulation
            ensemble.ensemble_grow(days=300)

            # Write ensemble results to output table or file
            ensemble.ensemble_results_to_output(database=metadata,
                                                pad_results=datetime.date(year,10,31))

            # Set status of current task to 'Finished'
            taskmanager.set_task_finished(task)

        except SQLAlchemyError, inst:
            print ("Database error: %s" % inst)
            # Break because of error in the database connection
            break

        except (dbi.SitedataError, dbi.SoildataError, dbi.TimerdataError,
                dbi.MeteodataError), inst:
            print "Error in WOFOST input procedures: %s" % inst
            # Set status of current task to 'Error'
            taskmanager.set_task_error(task)

        except Exception, inst:
            print ("General error in WOFOST procedures for year %i " + \
                     "crop %i and grid %i: %s") % (year, crop_no, grid_no, inst)
            # Set status of current task to 'Error'
            taskmanager.set_task_error(task)

        finally:
            #Get new task
            task = taskmanager.get_task()
    
    # Close files, logging and database connection
    logging.shutdown()
    connection.close()


def run_example_enkf(dsn=None):
    """Example for running PyWOFOST in ensemble Kalman filter (EnKF) mode
    
    Runs PyWOFOST for 6 crop types for one location in water-limited mode
    with the EnKF enambled using an ensemble of 50 members. Each member is
    initialized with different values for the initial soil moisture and each
    member receives different values for rainfall as forcing variable. External
    soil moisture estimates are retrieved from the table 'data_for_assimilation'
    and assimilated during the model run.     
    Executing PyWOFOST runs is done through the task manager. PyWOFOST output
    is writted to the database table 'pywofost_output', EnKF results like
    innovations are written to the table 'enkf_output'.
    
    Parameters:
    dsn - SQLAlchemy data source name pointing to the database to be used.
    """

    from .pywofost_enkf import PyWofostEnKF
    from . import db_input as dbi
    from .taskmanager import TaskManager
    from sqlalchemy.exceptions import SQLAlchemyError
    
    installdir = os.path.dirname(os.path.abspath(__file__))
    logconfig = os.path.join(installdir, 'logging.conf') 
    logging.config.fileConfig(logconfig)

    if dsn is None: # Assume SQlite demo DB
        db_location = os.path.join(installdir, "database","pywofost.db")
        dsn = "sqlite:///" + db_location

    # Open database connection and empty output table
    engine = create_engine(dsn)
    connection = engine.connect()
    metadata = MetaData(engine)
    table_pywofost_output = Table('pywofost_output', metadata,
                                  autoload=True)
    table_pywofost_output.delete().execute()
    table_enkf_output = Table('enkf_output', metadata, autoload=True)
    table_enkf_output.delete().execute()
    table_tasklist = Table("tasklist", metadata, autoload=True)
    table_tasklist.update().execute(status='Pending', hostname='')
    
    # Define ensemble size
    ensemble_size = 50

    # Initialise task manager
    taskmanager = TaskManager(metadata, connection, dbtype="SQLite",
                              hostname=socket.gethostname())
    # Loop until no tasks are left
    task = taskmanager.get_task()
    logger = None
    while task is not None:
        try:
            print "Running task: %i" % (task["task_id"])
            task_id = task["task_id"]
            grid_no = task["grid_no"]
            year = task["year"]
            crop_no = task["crop_no"]
            mode = task["sim_mode"]

            # Define random variable, take seed from database
            rv = random.RandomState(seed=int(task["randomseed"]))
        
            # Get Timer settings, Site, soil and Crop parameters
            timerdata = dbi.fetch_timerdata(metadata, grid_no, year, crop_no)
            sitedata  = dbi.fetch_sitedata(metadata, grid_no, year)
            cropdata = dbi.fetch_cropdata(metadata, grid_no, year, crop_no)
            soildata = dbi.fetch_soiltype(metadata, grid_no)
        
            # Get Meteo input data
            meteo_start_date = timerdata["WB_START_DATE"]
            meteo_end_date = timerdata["CROP_END_DATE"]
            meteof = dbi.MeteoFetcher(metadata, grid_no=grid_no,
                                     startdate=meteo_start_date,
                                     enddate=meteo_end_date,
                                     ensemble_mode=True)
                                     
            # Get data that should be assimilated.
            obs_params = {'soil_water_index':soildata}
            data_for_assimilation = dbi.fetch_assimdata(metadata, grid_no,
                                                crop_no, meteo_start_date,
                                                meteo_end_date, obs_params)

            # Generate ensemble of initial soil moisture values
            Esitedata = initial_soil_moisture(ensemble_size, rv, 
                                              soildata["SMW"],
                                              soildata["SMFCF"],
                                              soildata["SM0"],
                                              soildata["RDMSOL"])

            # Initialize ensemble
            ensemble = PyWofostEnKF(sitedata, timerdata, soildata, cropdata,
                                    meteof, data_for_assimilation, rv,
                                    ensemble_size=ensemble_size, mode=mode,
                                    Esitedata=Esitedata)
            ensemble.grow_with_assimilation()

            # Write ensemble results to output table
            ensemble.ensemble_results_to_output(database=metadata)
            print "Simulation results written to table 'pywofost_output'."

            ensemble.EnKF_results_to_output_device(database=metadata)
            print "EnKF results written to table 'enkf_output'."

            # Set status of current task to 'Finished'
            taskmanager.set_task_finished(task)

        except SQLAlchemyError, inst:
            print ("Database error: %s" % inst)
            # Break because of error in the database connection
            break

        except (dbi.SitedataError, dbi.SoildataError, dbi.TimerdataError,
            dbi.MeteodataError, dbi.DfaError), inst:
            print "Error in WOFOST input procedures: %s" % inst
             # Set status of current task to 'Error'
            taskmanager.set_task_error(task)

        except Exception, inst:
            print ("General error in WOFOST procedures for year %i " + \
                     "crop %i and grid %i: %s") % (year, crop_no, grid_no, inst)
            # Set status of current task to 'Error'
            taskmanager.set_task_error(task)

        finally:
            #Get new task
            task = taskmanager.get_task()
    
    # Close files, logging and database connection
    logging.shutdown()
    connection.close()

