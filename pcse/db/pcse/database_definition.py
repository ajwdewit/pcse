# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
################################
# WARNING: this code is obsolete. Left here as it has some useful table definitions.
################################
"""Create a PCSE database on any SQLAlchemy supported database

Creates the PCSE demo database structure on the target database and fills it
with a demo set of records from the SQLite database 'PCSE.db' that is
included with the PCSE distribution.

Public functions:
  migrate_db : migrates the PCSE database
  """
from __future__ import print_function
import os
from sqlalchemy import *

#-------------------------------------------------------------------------------
# Redefine SQLAlchemy Float() type because the standard version doesn't allow
# specifying the number of digits after the decimal point
class NewFloat(Numeric):
    def __init__(self, precision = 10, scale = 5, asdecimal=True):
        self.precision = precision
        self.scale = scale
    def adapt(self, impltype):
        return impltype(precision=self.precision, scale=self.scale)

#-------------------------------------------------------------------------------
def _create_PCSE_tables(metadata):
    """Creates the PCSE DB tables, needs an SQLAlchemy metadata object.
    """

    table_crop = Table('crop',metadata,
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('crop_name',String(length=40),nullable=False),
        Column('cropgroup_no',Integer,nullable=False),
        Column('crop_model',Integer,nullable=False),schema=None)
    table_crop.create()
    
    table_crop_calendar = Table('crop_calendar',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('variety_no',Integer,nullable=False),
        Column('year',Integer,primary_key=True,nullable=False),
        Column('start_type',String(length=20),nullable=False),
        Column('start_day',Date,nullable=False),
        Column('end_type',String(length=20),nullable=False),
        Column('end_day',Date,nullable=False),
        Column('max_duration',Integer,nullable=False),
        schema=None)
    table_crop_calendar.create()
    
    table_crop_parameter_value = Table('crop_parameter_value',metadata,
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('parameter_code',String(length=20),primary_key=True,nullable=False),
        Column('parameter_xvalue',NewFloat(),nullable=False),
        Column('parameter_yvalue',NewFloat()),schema=None)
    table_crop_parameter_value.create()
    
    table_grid = Table('grid',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('latitude',NewFloat(),nullable=False),
        Column('longitude',NewFloat(),nullable=False),
        Column('altitude',NewFloat(),nullable=False),
        Column('climate_barrier_no',Integer,nullable=False),
        Column('distance_to_coast',NewFloat(),nullable=False),schema=None)
    table_grid.create()
    
    table_grid_weather = Table('grid_weather',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('day',Date,primary_key=True,nullable=False),
        Column('maximum_temperature',NewFloat(),nullable=False),
        Column('minimum_temperature',NewFloat(),nullable=False),
        Column('vapour_pressure',NewFloat(),nullable=False),
        Column('windspeed',NewFloat(),nullable=False),
        Column('rainfall',NewFloat(),nullable=False),
        Column('e0',NewFloat(),nullable=False),
        Column('es0',NewFloat(),nullable=False),
        Column('et0',NewFloat(),nullable=False),
        Column('calculated_radiation',Integer(length=6),nullable=False),
        Column('snow_depth',NewFloat()),schema=None)
    table_grid_weather.create()
    
    table_pywofost_output = Table('pywofost_output',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('year',Integer,primary_key=True,nullable=False),
        Column('day',Date,primary_key=True,nullable=False),
        Column('simulation_mode',String(length=3),primary_key=True,nullable=False),
        Column('member_id',Integer,primary_key=True,nullable=False),
        Column('dvs',NewFloat(),nullable=True),
        Column('lai',NewFloat(),nullable=True),
        Column('tagp',NewFloat(),nullable=True),
        Column('twso',NewFloat(),nullable=True),
        Column('twlv',NewFloat(),nullable=True),
        Column('twst',NewFloat(),nullable=True),
        Column('twrt',NewFloat(),nullable=True),
        Column('sm',NewFloat(),nullable=True),
        Column('tra',NewFloat(),nullable=True),
        Column('rd',NewFloat(),nullable=True),
        schema=None)
    table_pywofost_output.create()

    table_pywofost_unittest = Table('pywofost_unittest_benchmarks',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('year',Integer,primary_key=True,nullable=False),
        Column('day',Date,primary_key=True,nullable=False),
        Column('simulation_mode',String(length=3),primary_key=True,nullable=False),
        Column('member_id',Integer,primary_key=True,nullable=False),
        Column('dvs',NewFloat(),nullable=False),
        Column('lai',NewFloat(),nullable=False),
        Column('tagp',NewFloat(),nullable=False),
        Column('twso',NewFloat(),nullable=False),
        Column('twlv',NewFloat(),nullable=False),
        Column('twst',NewFloat(),nullable=False),
        Column('twrt',NewFloat(),nullable=False),
        Column('sm',NewFloat(),nullable=False),
        Column('tra',NewFloat(),nullable=False),
        Column('rd',NewFloat(),nullable=False),
        schema=None)
    table_pywofost_unittest.create()
    
    table_enkf_output = Table('enkf_output', metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('year',Integer,primary_key=True,nullable=False),
        Column('day',Date,primary_key=True,nullable=False),
        Column('variable',String(length=255),nullable=False),
        Column('A', String(length=1000), nullable=False),
        Column('D', String(length=1000), nullable=False),
        Column('P_e', String(length=1000), nullable=False),
        Column('R_e', String(length=1000), nullable=False),
        Column('K', String(length=1000), nullable=False),
        schema=None)
    table_enkf_output.create()
    
    table_rooting_depth = Table('rooting_depth',metadata,
        Column('rooting_depth_class',Integer,nullable=False),
        Column('maximum_rootable_depth',Integer,nullable=False),
        schema=None)
    table_rooting_depth.create()
    
    table_site = Table('site',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('year',Integer,primary_key=True,nullable=False),
        Column('ifunrn',Integer,nullable=False),
        Column('max_surface_storage',NewFloat(),nullable=False),
        Column('not_infiltrating_fraction',NewFloat(),nullable=False),
        Column('initial_surface_storage',NewFloat(),nullable=False),
        Column('inital_water_availability',NewFloat(),nullable=False),
        schema=None)
    table_site.create()
    
    table_soil_physical_group = Table('soil_physical_group',metadata,
    
        Column('soil_group_no',Integer,primary_key=True,nullable=False),
        Column('parameter_code',String(length=30),primary_key=True,nullable=False),
        Column('parameter_xvalue',NewFloat(),nullable=False),
        Column('parameter_yvalue',NewFloat(),nullable=True),
        schema=None)
    table_soil_physical_group.create()
    
    table_soil_type = Table('soil_type',metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('rooting_depth_class',Integer,primary_key=True,nullable=False),
        Column('soil_group_no',Integer,primary_key=True,nullable=False),
        schema=None)
    table_soil_type.create()
    
    table_variety_parameter_value = Table('variety_parameter_value',metadata,
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('variety_no',Integer,primary_key=True,nullable=False),
        Column('parameter_code',String(length=20),primary_key=True,nullable=False),
        Column('parameter_xvalue',NewFloat(),nullable=False),
        Column('parameter_yvalue',NewFloat()),schema=None)
    table_variety_parameter_value.create()

    table_data_for_assimilation = Table('data_for_assimilation', metadata,
        Column('grid_no',Integer,primary_key=True,nullable=False),
        Column('crop_no',Integer,primary_key=True,nullable=False),
        Column('day',Date,primary_key=True,nullable=False),
        Column('observed_state',String(length=25),primary_key=True,nullable=False),
        Column('value',NewFloat(),nullable=False),
        Column('variance',NewFloat(),nullable=False),
        schema=None)
    table_data_for_assimilation.create()
    
    table_egw = Table('ensemble_grid_weather',metadata,
        Column('grid_no',Integer,primary_key=True),
        Column('day',Date,primary_key=True),
        Column('ensemble_id',Integer,primary_key=True),
        Column('rainfall',NewFloat(),nullable=False),
        schema=None)
    table_egw.create()
    
    table_tasklist = Table('tasklist',metadata,
        Column('task_id',Integer,primary_key=True, nullable=False),
        Column('status',String(length=15),nullable=False),
        Column('hostname',String(length=50),nullable=True),
        Column('grid_no',Integer, nullable=False),
        Column('crop_no',Integer,nullable=False),
        Column('year',Integer,nullable=False),
        Column('sim_mode',String(length=5),nullable=True),
        Column('randomseed',Integer,nullable=False),
        schema=None)
    table_tasklist.create()
        
    print("Tables created!")

#-------------------------------------------------------------------------------
def _fill_pywofost_tables(engine, metadata, table_collection):
    """Fills the pywofost tables with data.
    
    Needs an SQLAlchemy metadata object and the table_collection
    (e.g. {tablename:records, ...}) whose data is to be inserted into the
    new PyWofost database.
    """
    
    exc = []
    for tablename in table_collection:
        print("Filling table %s" % tablename)
        if len(table_collection[tablename]) > 0:
            t = Table(tablename, metadata, autoload=True)
            try:
                engine.execute(t.insert(), table_collection[tablename])
            except Exception as inst:
                exc += [inst]
    for e in exc:
        print(e)

#-------------------------------------------------------------------------------
def _retrieve_records_from_source(metadata):
    """Retrieves all records from all input tables in source DB"""
    
    input_tables = ['crop', 'crop_calendar', 'grid',
                'grid_weather', 'pywofost_unittest_benchmarks',
                'rooting_depth', 'soil_physical_group', 'soil_type',
                'variety_parameter_value','data_for_assimilation',
                'ensemble_grid_weather','tasklist', 'site',
                'crop_parameter_value']
    table_collection = {}
    for table in input_tables:
        t = Table(table, metadata, autoload=True)
        r = select([t]).execute()
        records = []
        for rec in r:
            records += [dict(rec)]
        table_collection[table] = records
        
    return table_collection
        
#-------------------------------------------------------------------------------
def migrate_db(target_dsn=None):
    """Migrates the structure and data in the PyWOFOST demo database.

keyword parameters: 
    target_dsn : SQLAlchemy connection string specifying the database to 
                 connect to.

Examples of SQLalchemy connection strings"
    For MySQL
      target_dsn = "mysql://<user>:<password>@<hostname>/<database>"
    For ORACLE
      target_dsn = "oracle://<user>:<password>@TNS"
    For SQLite (on windows)
      target_dsn = "sqlite:///D:/DATA/pywofost.db"
    For SQLite (on UNIX)
      target_dsn = "sqlite:////home/user/pywofost.db"

See: http://www.sqlalchemy.org/docs/06/core/engines.html#supported-databases
     for more examples and other supported databases
"""

    # Open target database connection
    if target_dsn is None:
        print("No target_dsn specified, see docstring on migrate_db() function"
              " for dsn specification.")
        return
    try:
        target_engine = create_engine(target_dsn)
        target_metadata = MetaData(target_engine)
    except Exception as e:
        print("Unable to open connection to database, due to following exception:"
              " %s" % e.args[0])
        return

    # Open source database connection
    installdir = os.path.dirname(os.path.abspath(__file__))
    db_location = os.path.join(installdir, "pywofost.db")
    dsn = "sqlite:///" + db_location
    try:
        source_engine = create_engine(dsn)
        source_metadata = MetaData(source_engine)
    except Exception as e:
        print("Unable to open connection to pywofost demo database"
              "with the following exception:\n %s" % e.args[0])
    
    # Build database
    try:
        _create_pywofost_tables(target_metadata)
        table_collection = _retrieve_records_from_source(source_metadata)
        _fill_pywofost_tables(target_engine, target_metadata, table_collection)
        print("PyWOFOST demo database succesfully migrated!")
    except Exception as e:
        print("Migrating the database failed with the "
              "following exception:\n %s" % e.args[0])

