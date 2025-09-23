# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), March 2024
"""Module defines code for running unittests on the complete PCSE/WOFOST model.

Unit tests for individual components are included in the sources files of these
components.

Classes defined here:
* WofostBenchmarkRetriever
* WofostOutputRetriever
* WofostTestingTemplate (and derived classes for testing)

Functions defined here:
* run_units_tests

"""
import os
import random
import unittest

import pandas as pd
import sqlite3

from .run_wofost import run_wofost
from ..settings import settings


def dict_factory(cursor, row):
    fields = [column[0] for column in cursor.description]
    return {key: value for key, value in zip(fields, row)}


class WofostBenchmarkRetriever:
    """Retrieves benchmark results for PCSE WOFOST.
    
    Class retrieves benchmarks from the PCSE DB that are being used to
    unittest (benchmark) the WOFOST output. The benchmarks are stored in the
    table 'wofost_unittest_benchmarks'.
    
    Example:
    retriever = WofostBenchmarkRetriever('data source name', crop, grid, mode)
    benchmark_data = retriever('development_stage')
    """

    def __init__(self, dsn, crop, grid, mode='pp'):
        DBconn = sqlite3.connect(dsn)
        DBconn.row_factory = dict_factory
        sql = f"select * from wofost_unittest_benchmarks where crop_no=? and grid_no=? and " \
              f"simulation_mode=? and member_id=?"
        cursor = DBconn.cursor()
        r = cursor.execute(sql, (crop, grid, mode, 0))
        df = pd.DataFrame(r.fetchall())
        self.df_benchmarks = df.set_index("day")
        DBconn.close()
        self.crop = crop
        self.grid = grid
        self.mode = mode
        self.member_id = 0
    
    def __call__(self, variable):
        """Retrieves benchmark data for specified variable: [(day, variable),..]
        """
        if not variable in self.df_benchmarks.columns:
            msg = f"variable {variable} not available in DataFrame!"
            raise RuntimeError(msg)

        rows = self.df_benchmarks[[variable]].to_records()
        return rows


class WofostOutputRetriever:
    """Retrieves results from a Wofost simulation.
    
    Class retrieves results from a Wofost simulation based on the table
    'sim_results_timeseries'. These results are then compared to the benchmark
    results for unit testing. This procedure assumes that there is only a single
    simulation (grid, crop, year) present in the table 'sim_results_timeseries'
    as no selection on (grid, crop, year) is being done.
    
    Example:
    retriever = WofostOutputRetriever('data source name')
    
    one_day = retriever(date(2000,1,1), 'development_stage')
    last_day = retriever.getWofostOutputLastDay('development_stage')
    """

    def __init__(self, dsn):
        DBconn = sqlite3.connect(dsn)
        DBconn.row_factory = dict_factory
        sql = f"select * from sim_results_timeseries"
        cursor = DBconn.cursor()
        r = cursor.execute(sql)
        df = pd.DataFrame(r.fetchall())
        self.df_sim_results = df.set_index("day")
        self.maxday = self.df_sim_results.index.max()
        DBconn.close()

    def __call__(self, day, variable):
        """Returns the specified WOFOST variable for specified day."""

        ix = self.df_sim_results.index == day
        if not ix.any():
            msg = f"cannot find simulation results for day {day} and variable {variable}"
            raise RuntimeError(msg)
        value = self.df_sim_results.loc[ix, variable][0]
        return float(value)


class WofostTestingTemplate(unittest.TestCase):
    """Template for executing WOFOST unit tests.
    
    The template defines the setUp() and runTest() routines that are common to
    all WOFOST unit testing runs. Most of functionality comes from
    subclassing 'unittest.TestCase'. Note that each unittest is simply defined
    as a subclass of this template. Only the crop, grid, year and mode are
    test-specific and thus defined as test class attributes.
    
    To prevent false FAILS caused by different python versions, databases and
    cpu architectures, biomass values are only checked for accuracy up to one
    or three decimals.    
    """

    benchmark_vars = [("DVS",3), ("TRA",3), ("RD",3),("SM", 3), ("LAI",2),
                      ("TAGP",1),("TWLV",1),("TWST",1),("TWSO",1),("TWRT",1)]

    def __init__(self, testname):
        db_location = os.path.join(settings.PCSE_USER_HOME, "pcse.db")
        self.dsn = os.path.normpath(db_location)
        unittest.TestCase.__init__(self, testname)
        
    def setUp(self):
        "Sets up the simulation in order to verify the results"
        run_wofost(dsn=self.dsn, crop=self.crop, year=self.year,
                     grid=self.grid, mode=self.mode, clear_table=True)
        self.OutputRetriever = WofostOutputRetriever(dsn=self.dsn)
        self.BenchmarkRetriever = WofostBenchmarkRetriever(dsn=self.dsn,
                                                             crop=self.crop,
                                                             grid=self.grid,
                                                             mode=self.mode)

    def run_benchmark(self, var_to_benchmark, precision):
        benchmark_data = self.BenchmarkRetriever(var_to_benchmark)
        msg = "Failure to retrieve benchmark data."
        self.assertTrue(len(benchmark_data) > 0, msg)
        n_assert = 0
        for (day, benchmark) in benchmark_data:
            value = self.OutputRetriever(day, var_to_benchmark)
            if value is None:
                continue
            diff = float(abs(value - benchmark))
            assertmsg = "Test day, variable %s, %s: %f vs %f" % \
                        (day, var_to_benchmark, benchmark, value)
            self.assertAlmostEqual(diff, 0., precision, assertmsg)
            n_assert += 1
        # Test if any assertions were done at all which can occur
        # due to missing records (e.g dekadal or monthly output),
        # a mis-spelled variable name, or a full series of records with
        # NULL (e.g. None) values in the database.
        msg = "No data found in sim_results_timeseries table for variable '%s'"
        msg = msg % var_to_benchmark
        self.assertGreater(n_assert, 0, msg)

    def runTest(self):
        for eachvar, precision in self.benchmark_vars:
            self.run_benchmark(eachvar, precision)


class TestPotentialWinterWheat(WofostTestingTemplate):
    crop = 1
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedWinterWheat(WofostTestingTemplate):
    crop = 1
    year = 2000
    grid = 31031
    mode = "wlp"


class TestPotentialGrainMaize(WofostTestingTemplate):
    crop = 2
    year = 2000
    grid = 31031
    mode = "pp"


class TestWaterlimitedGrainMaize(WofostTestingTemplate):
    crop = 2
    year = 2000
    grid = 31031
    mode = "wlp"


class TestPotentialSpringBarley(WofostTestingTemplate):
    crop = 3
    year = 2000
    grid = 31031
    mode = "pp"


class TestWaterlimitedSpringBarley(WofostTestingTemplate):
    crop = 3
    year = 2000
    grid = 31031
    mode = "wlp"


class TestPotentialPotato(WofostTestingTemplate):
    crop = 7
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedPotato(WofostTestingTemplate):
    crop = 7
    year = 2000
    grid = 31031
    mode = "wlp"


class TestPotentialWinterRapeseed(WofostTestingTemplate):
    crop = 10
    year = 2000
    grid = 31031
    mode = "pp"


class TestWaterlimitedWinterRapeseed(WofostTestingTemplate):
    crop = 10
    year = 2000
    grid = 31031
    mode = "wlp"


class TestPotentialSunflower(WofostTestingTemplate):
    crop = 11
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedSunflower(WofostTestingTemplate):
    crop = 11
    year = 2000
    grid = 31031
    mode = "wlp"


def suite():
    """returns the unit tests for the PCSE/WOFOST model.
    
    keyword parameters:
    dsn : SQLAlchemy data source name for the database that should be used.
    """

    suite = unittest.TestSuite()
    tests = [TestPotentialWinterWheat('runTest'),
             TestWaterlimitedWinterWheat('runTest'),
             TestPotentialGrainMaize('runTest'),
             TestWaterlimitedGrainMaize('runTest'),
             TestPotentialSpringBarley('runTest'),
             TestWaterlimitedSpringBarley('runTest'),
             TestPotentialPotato('runTest'),
             TestWaterlimitedPotato('runTest'),
             TestPotentialWinterRapeseed('runTest'),
             TestWaterlimitedWinterRapeseed('runTest'),
             TestPotentialSunflower('runTest'),
             TestWaterlimitedSunflower('runTest')]

    # Shuffle tests in random order in order to test whether different results
    # are obtained when tests are run in a different order.
    random.shuffle(tests)
    suite.addTests(tests)
    return suite
