# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
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

from sqlalchemy import create_engine, MetaData, select, Table, and_, func
from sqlalchemy.types import Date as saDate

from .run_wofost import run_wofost
from ..settings import settings


#-------------------------------------------------------------------------------
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
        engine = create_engine(dsn)
        meta = MetaData(engine)
        self.table_benchm = Table('wofost_unittest_benchmarks',
                                  meta, autoload=True)
        self.crop = crop
        self.grid = grid
        self.mode = mode
        self.member_id = 0
    
    def __call__(self, variable):
        "Retrieves benchmark data for specified variable: [(day, variable),..]"
        t1 = self.table_benchm.alias()
        column_for_retrieval = t1.c[variable]
        if column_for_retrieval is not None:
            s = select([t1.c.day, column_for_retrieval],
                       and_(t1.c.crop_no==self.crop,
                            t1.c.grid_no==self.grid,
                            t1.c.simulation_mode==self.mode,
                            t1.c.member_id==self.member_id))
            rows = s.execute().fetchall()
            return rows
        else:
            raise KeyError("Variable not present in sim_results_timeseries table.")

#-------------------------------------------------------------------------------
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
        db_engine = create_engine(dsn)
        db_metadata = MetaData(db_engine)
        self.table_pwo = Table('sim_results_timeseries', db_metadata,
                               autoload=True)
        s = select([func.max(self.table_pwo.c.day, type_=saDate)])
        self.maxday = s.execute().fetchone()[0]
    
    def __call__(self, day, variable):
        """Returns the specified WOFOST variable for specified day."""
        s = select([self.table_pwo], and_(self.table_pwo.c.day==day))
        row = s.execute().fetchone()
        #print "day %s" % day
        return row[variable]
        
    def getWofostOutputLastDay(self, variable):
        """Returns the specified WOFOST variable on the last day."""
        s = select([self.table_pwo], and_(self.table_pwo.c.day==self.maxday))
        row = s.execute().fetchone()
        return row[variable]

#-------------------------------------------------------------------------------
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

    def __init__(self, testname, dsn=None):
        if dsn is None: # Assume SQlite demo DB
            db_location = os.path.join(settings.PCSE_USER_HOME,"pcse.db")
            db_location = os.path.normpath(db_location)
            dsn = "sqlite:///" + db_location
        self.dsn = dsn
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

class TestWaterlimitedSpringBarley1965(WofostTestingTemplate):
    crop = 3
    year = 1965     # CAMPAIGNYEAR
    grid = 1        # iso 35042
    mode = "wlp"
    multi_layer = True

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

def suite(dsn=None):
    """returns the unit tests for the PCSE/WOFOST model.
    
    keyword parameters:
    dsn : SQLAlchemy data source name for the database that should be used.
    """

    suite = unittest.TestSuite()
    tests = [TestPotentialWinterWheat('runTest', dsn),
             TestWaterlimitedWinterWheat('runTest', dsn),
             TestPotentialGrainMaize('runTest', dsn),
             TestWaterlimitedGrainMaize('runTest', dsn),
             TestPotentialSpringBarley('runTest', dsn),
             TestWaterlimitedSpringBarley('runTest', dsn),
             TestPotentialPotato('runTest', dsn),
             TestWaterlimitedPotato('runTest', dsn),
             TestPotentialWinterRapeseed('runTest', dsn),
             TestWaterlimitedWinterRapeseed('runTest', dsn),
             TestPotentialSunflower('runTest', dsn),
             TestWaterlimitedSunflower('runTest', dsn)]
    # tests = [TestPotentialWinterWheat('runTest', dsn),
    #         TestWaterlimitedWinterWheat('runTest', dsn),
    #         TestPotentialGrainMaize('runTest', dsn),
    #         TestPotentialSpringBarley('runTest', dsn),
    #         TestPotentialPotato('runTest', dsn),
    #         TestPotentialWinterRapeseed('runTest', dsn),
    #         TestPotentialSunflower('runTest', dsn)]
    # tests = [TestWaterlimitedSpringBarley1965('runTest', dsn)]

    # Shuffle tests in random order in order to test whether different results
    # are obtained when tests are run in a different order.
    random.shuffle(tests)
    suite.addTests(tests)
    return suite
