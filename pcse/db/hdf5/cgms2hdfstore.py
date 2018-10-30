# -*- coding: utf-8 -*-
# Copyright (c) 2004-2017 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), November 2017
"""
Converts the tables for running WOFOST from a CGMS12 database to a
pandas HDF5 store. This provides a very high performance for data retrieval.

See also: https://pandas.pydata.org/pandas-docs/stable/io.html#io-hdf5
"""
from __future__ import print_function
import argparse
import warnings
import sys, os
import time

import sqlalchemy as sa
from sqlalchemy import Table, and_, MetaData, create_engine, select
import pandas as pd
import numpy as np
import tables

from pcse.util import wind10to2


def create_parser():
    parser = argparse.ArgumentParser(description='Convert part of CGMS12 DB structure to HDF5.')
    parser.add_argument('--h5f', dest='h5fname', action='store', default=None,
                        help='Name of HDF5 Container to write to.')
    parser.add_argument('--dsn', dest='dsn', action='store', default=None,
                        help="SQLAlchemy connection URL for CGMS DB to connect to. See also "
                             "http://docs.sqlalchemy.org/en/latest/core/engines.html"
                        )
    parser.add_argument('--include_grid_weather', dest='incl_gw', action='store_true',
                        help="Include the table with gridded weather in the HDF5 file."
                        )
    parser.add_argument('--include_initial_soil_water', dest='incl_isw', action='store_true',
                        help="Include the table with values for initial soil water in the HDF5 file."
                        )
    parser.add_argument('--include_crop_soil_tables', dest='incl_cropsoil', action='store_true',
                        help="Include the tables related to crop and soil parameters in the HDF5 file."
                        )
    parser.add_argument('--include_crop_calendar', dest='incl_calendar', action='store_true',
                        help="Include the crop calendar table in the HDF5 file."
                        )
    parser.add_argument('--include_all', dest='incl_all', action='store_true',
                        help="Include all tables in the HDF5 file (file size may become very large!)."
                        )
    parser.add_argument('--crop_no', dest='crop_no', action='store', default=None, type=int,
                        help="Only include records for given crop_no. This can reduce file size considerably."
                        )
    return parser


def fetch_grids(engine):
    """retrieves the GRID table from a CGMS12 database.
    """
    meta = MetaData(engine)
    tbl_grid = Table("grid", meta, autoload=True)
    s = tbl_grid.select()
    df_grids = pd.read_sql(s, engine)
    return df_grids


def store_small_tables(engine, fname_HDFstore):
    """retrieves the a set of smaller tables from a CGMS12 database.
    """
    meta = MetaData(engine)
    tables = {"suitability": ["cropgroup_no"],
              "soil_typologic_unit": ["stu_no"],
              "rooting_depth": ["class"],
              "soil_physical_group": ["soil_group_no"],
              "emu": ["grid_no","smu_no"],
              "soil_association_composition":["smu_no", "stu_no"],
              "crop": ["crop_no"],
              "crop_parameter_value":["crop_no"],
              "variety_parameter_value":["crop_no", "variety_no"],
              "site":[]
              }

    with pd.io.pytables.HDFStore(fname_HDFstore) as store:
        for tbl_name, dcolumns in tables.items():
            print("Storing table: %s" % tbl_name)
            tbl = Table(tbl_name, meta, autoload=True)
            df = pd.read_sql(tbl.select(), engine)
            store.put("/%s" % tbl_name, df, format="table", data_columns=dcolumns)


def store_crop_calendar(engine, fname_HDFstore, args):
    """retrieves the CROP_CALENDAR table from a CGMS12 database.

    if the --crop_no option is used only the records for the given crop_no
    will be retrieved.
    """
    meta = MetaData(engine)
    tbl_cal = Table("crop_calendar", meta, autoload=True)

    # retrieve distinct crop types from DB table
    s = sa.select([tbl_cal.c.crop_no]).distinct()
    crops = [row[0] for row in s.execute()]

    # Check if only one crop type should be selected from DB
    if args.crop_no is not None:
        if args.crop_no not in crops:
            print("Crop ID specified with --cropno (%s) not found in CROP_CALENDAR table! Returning..." % args.crop_no)
            sys.exit()
        crops = [args.crop_no]

    # Start pulling crop_calendar data from DB
    dataset_name = "/crop_calendar"
    with pd.io.pytables.HDFStore(fname_HDFstore) as store:
        for crop in crops:
            print("Storing crop_calendar for crop %i" % crop)
            s = tbl_cal.select().where(tbl_cal.c.crop_no==crop)
            df_cal = pd.read_sql(s, engine)
            if dataset_name in store:
                store.append(dataset_name, df_cal, data_columns=["grid_no", "crop_no", "year"])
            else:
                store.put(dataset_name, df_cal, format="table", data_columns=["grid_no", "crop_no", "year"])


def store_initial_soil_water(engine, fname_HDFstore, args):
    """retrieves the INITIAL_SOIL_WATER table from a CGMS12 database.

    if the --crop_no option is used only the records for the given crop_no
    will be retrieved.
    """
    meta = MetaData(engine)
    tbl_isw = Table("initial_soil_water", meta, autoload=True)

    # retrieve distinct crop types from DB table
    if args.crop_no is not None:
        s = sa.select([tbl_isw.c.crop_no]).distinct()
        crops = [row[0] for row in s.execute()]
        if args.crop_no not in crops:
            print("Crop ID specified with --cropno (%s) not found in INITIAL_SOIL_WATER table! Returning..." % args.crop_no)
            sys.exit()

    # Select distinct years to iterate over
    s = sa.select([tbl_isw.c.year]).distinct()
    years = s.execute()
    dataset_name = "/initial_soil_water"
    with pd.io.pytables.HDFStore(fname_HDFstore) as store:
        for yr, in sorted(years):
            if args.crop_no:
                s = tbl_isw.select().where(sa.and_(tbl_isw.c.year == yr,
                                                   tbl_isw.c.crop_no == args.crop_no))
                print("Storing initial_soil_water for crop %i and year %i" % (args.crop_no, yr))
            else:
                s = tbl_isw.select().where(tbl_isw.c.year == yr)
                print("Storing initial_soil_water for year %i" % yr)
            df_isw = pd.read_sql(s, engine)
            if dataset_name in store:
                store.append(dataset_name, df_isw, data_columns=["grid_no", "stu_no", "crop_no", "year"])
            else:
                store.put(dataset_name, df_isw, format="table", data_columns=["grid_no", "stu_no", "crop_no", "year"])


def store_gridded_weather(engine, fname_HDFstore, args):
    """retrieves the WEATHER_OBS_GRID table from a CGMS12 database.

    The meteo data are stored under their own group to speed up the data
    retrieval and can be found under /<grid_no>/data
    """
    df_grids = fetch_grids(engine)
    ngrids = len(df_grids)
    warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)
    with pd.io.pytables.HDFStore(fname_HDFstore) as store:
        store.put('/grid', df_grids, format='table', data_columns=True)
        for grid_row in df_grids.itertuples():
            grid_done = False
            while not grid_done:
                try:
                    if engine is None:
                        engine = connect_to_db(args.dsn)

                    if grid_row.Index % 100 == 0:
                        print("\nProcessing %i out of %i grids" % (grid_row.Index, ngrids))
                    else:
                        print(".", end="")
                        sys.stdout.flush()  # directly print "." to terminal
                    df_meteo_raw = fetch_weather_from_db(engine, int(grid_row.grid_no))
                    df_meteo_pro = process_meteo(df_meteo_raw)
                    store.put('%i/data' % grid_row.grid_no, df_meteo_pro)
                    grid_done = True
                except sa.exc.SQLAlchemyError as e:  # DB connection failure
                    engine = None
                    print("DB connection failure: %s" % e)
                except Exception as e:
                    print("General exception - aborting: %s" % e)
                    sys.exit()


def fetch_weather_from_db(engine, grid_no):
    """Retrieves the meteo data from table 'weather_obs_grid'
    """

    # if start_date/end_date are None, define a date in the far past/future
    meta = MetaData(engine)
    table_db = Table("weather_obs_grid", meta, autoload=True)
    s = select([table_db], and_(table_db.c.grid_no == grid_no)
               )
    df_meteo = pd.read_sql(s, engine)
    df_meteo = df_meteo.set_index("day")
    return df_meteo


def process_meteo(df_meteo):
    """ Converts the meteo table into a structure directly usable by PCSE.
    :param df_meteo: a data frame with records from a CGMS12 WEATHER_OBS_GRID
    :return: a new dataframe with the appropriate structure.
    """

    df_new = pd.DataFrame({"TMAX": df_meteo.temperature_max.astype(np.float32),
                           "TMIN": df_meteo.temperature_min.astype(np.float32),
                           "TEMP": df_meteo.temperature_avg.astype(np.float32),
                           "VAP": df_meteo.vapourpressure.astype(np.float32),
                           "WIND": df_meteo.windspeed.apply(wind10to2).astype(np.float32),
                           "RAIN": (df_meteo.precipitation / 10.).astype(np.float32),
                           "IRRAD": (df_meteo.radiation * 1000.).astype(np.float32),
                           "SNOWDEPTH": df_meteo.snowdepth.astype(np.float32),
                           "E0": (df_meteo.e0 / 10.).astype(np.float32),
                           "ES0": (df_meteo.es0 / 10.).astype(np.float32),
                           "ET0": (df_meteo.et0 / 10.).astype(np.float32)},
                          index=df_meteo.index)
    return df_new


def connect_to_db(dsn):
    retries = 0
    while True:
        try:
            engine = create_engine(dsn)
            engine.connect()
            break
        except Exception as e:
            if retries == 0:
                print("DB connection failed, retrying", end="")
            else:
                print(".", end="")
                sys.stdout.flush()
            retries += 1
            if retries > 5:
                print("\nDB connection failed after 5 retries: %s" % e)
                sys.exit()
            time.sleep(5)
    return engine


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not args.h5fname or not args.dsn:
        parser.print_help()
        return

    fname_HDFStore = os.path.abspath(args.h5fname)
    engine = connect_to_db(args.dsn)

    print("Writing data to HDF store: %s" % fname_HDFStore)
    print("Fetching data from DB: %s" % engine)

    done = []
    if args.incl_cropsoil or args.incl_all:
        store_small_tables(engine, fname_HDFStore)
        done.append("crop_soil")
    if args.incl_calendar or args.incl_all:
        store_crop_calendar(engine, fname_HDFStore, args)
        done.append("crop_calendar")
    if args.incl_isw or args.incl_all:
        store_initial_soil_water(engine, fname_HDFStore, args)
        done.append("initial_soil_water")
    if args.incl_gw or args.incl_all:
        store_gridded_weather(engine, fname_HDFStore, args)
        done.append("gridded weather")

    if not done:
        msg = """No tables included for conversion to HDF5, nothing done...
        See --include options!"""
        print(msg)
        parser.print_help()


if __name__ == "__main__":
    main()