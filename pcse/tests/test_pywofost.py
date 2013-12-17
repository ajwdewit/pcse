"""Module defines code for running unittests on the complete PyWofost model.

Unit tests for individual components are included in the sources files of these
components.

Classes defined here:
* PyWofostBenchmarkRetriever
* PyWofostOutputRetriever
* PyWofostTestingTemplate (and derived classes for testing)

Functions defined here:
* run_units_tests

"""
import os
import unittest
import random

from sqlalchemy import create_engine, MetaData, select, Table, and_, func
from sqlalchemy.types import Date as saDate

from ..run_pywofost import run_pywofost

#-------------------------------------------------------------------------------
class PyWofostBenchmarkRetriever:
    """Retrieves benchmark results for pywofost.
    
    Class retrieves benchmarks from the pywofost DB that are being used to
    unittest (benchmark) the PyWofost output. The benchmarks are stored in the
    table 'pywofost_unittest_benchmarks'.
    
    Example:
    retriever = PyWofostBenchmarkRetriever('data source name', crop, grid, mode)                
    benchmark_data = retriever('development_stage')
    """

    def __init__(self, dsn, crop, grid, mode='pp'):
        pywofost_engine = create_engine(dsn)
        pywofost_metadata = MetaData(pywofost_engine)
        self.table_benchm = Table('pywofost_unittest_benchmarks',
                                  pywofost_metadata, autoload=True)
        self.crop = crop
        self.grid = grid
        self.mode = mode
        self.member_id = 0
    
    def __call__(self, variable):
        "Retrieves benchmark data for specified variable: [(day, variable),..]"
        t1 = self.table_benchm.alias()
        column_for_retrieval = t1.c[variable]
        benchmark_data = []
        if column_for_retrieval is not None:
            s = select([t1.c.day, column_for_retrieval],
                       and_(t1.c.crop_no==self.crop,
                            t1.c.grid_no==self.grid,
                            t1.c.simulation_mode==self.mode,
                            t1.c.member_id==self.member_id))
            rows = s.execute().fetchall()
            return rows
        else:
            raise KeyError("Variable not present in pywofost_output table.")

#-------------------------------------------------------------------------------
class PyWofostOutputRetriever:
    """Retrieves results from a PyWofost simulation.
    
    Class retrieves results from a PyWofost simulation based on the table
    'pywofost_output'. These results are then compared to the benchmark results
    for unit testing. This procedure assumes that there is only a single
    simulation (grid, crop, year) present in the table 'pywofost_output'
    as no selection on (grid, crop, year) is being done.
    
    Example:
    retriever = PyWofostOutputRetriever('data source name')
    
    one_day = retriever(date(2000,1,1), 'development_stage')
    last_day = retriever.getPyWofostOutputLastDay('development_stage')    
    """

    def __init__(self, dsn):
        pywofost_engine = create_engine(dsn)
        pywofost_metadata = MetaData(pywofost_engine)
        self.table_pwo = Table('pywofost_output', pywofost_metadata,
                               autoload=True)
        s = select([func.max(self.table_pwo.c.day, type_=saDate)])
        self.maxday = s.execute().fetchone()[0]
    
    def __call__(self, day, variable):
        "Returns the specified pywofost variable for specified day."
        s = select([self.table_pwo], and_(self.table_pwo.c.day==day))
        row = s.execute().fetchone()
        #print "day %s" % day
        return row[(variable.lower())]
        
    def getPyWofostOutputLastDay(self, variable):
        "Returns the specified pywofost variable on the last day."
        s = select([self.table_pwo], and_(self.table_pwo.c.day==self.maxday))
        row = s.execute().fetchone()
        return row[(variable.lower())]

#-------------------------------------------------------------------------------
class PyWofostTestingTemplate(unittest.TestCase):
    """Template for executing PyWofost unit tests.
    
    The template defines the setUp() and runTest() routines that are common to
    all PyWofost unit testing runs. Most of functionality comes from
    subclassing 'unittest.TestCase'. Note that each unittest is simply defined
    as a subclass of this template. Only the crop, grid, year and mode are
    test-specific and thus defined as test class attributes.
    
    To prevent false FAILS caused by different python versions, databases and
    cpu architectures, biomass values are only checked for accuracy up to one
    or three decimals.    
    """

    benchmark_vars = [("dvs",3), ("tra",3), ("rd",3),("sm", 3), ("lai",2),
                      ("tagp",1),("twlv",1),("twst",1),("twso",1),("twrt",1)]
    multi_layer = False
    #benchmark_vars = [("dvs",3), ("lai",3), ("tra",3), ("rd",3),("sm", 3),
    #                  ("tagp",1),("twlv",1),("twst",1),("twso",1),("twrt",1)]
    # 2012/02/02: leave out rooting depth (rd) for benchmarking potential
    # production
    #benchmark_vars = [("dvs",3), ("lai",3), ("tra",3), ("sm", 3), 
    #                  ("tagp",1),("twlv",1),("twst",1),("twso",1),("twrt",1)]
    
    def __init__(self, testname, dsn=None):
        if dsn is None: # Assume SQlite demo DB
            installdir = os.path.dirname(os.path.abspath(__file__))
            db_location = os.path.join(installdir, "..", "db","pcse","pywofost.db")
            db_location = os.path.normpath(db_location)
            dsn = "sqlite:///" + db_location
        self.dsn=dsn
        unittest.TestCase.__init__(self, testname)
        
    def setUp(self):
        "Sets up the simulation in order to verify the results"
        run_pywofost(dsn=self.dsn, crop=self.crop, year=self.year,
                     grid=self.grid, mode=self.mode, clear_table=True,
                     multi_layer=self.multi_layer)
        self.OutputRetriever = PyWofostOutputRetriever(dsn=self.dsn)
        self.BenchmarkRetriever = PyWofostBenchmarkRetriever(dsn=self.dsn,
                                                             crop=self.crop,
                                                             grid=self.grid,
                                                             mode=self.mode)
    def run_benchmark(self, var_to_benchmark, precision):
        benchmark_data = self.BenchmarkRetriever(var_to_benchmark)
        msg = "Failure to retrieve benchmark data."
        self.failUnless(len(benchmark_data)>0, msg)
        for (day, benchmark) in benchmark_data:
            value = self.OutputRetriever(day, var_to_benchmark)
            if value is None:
                continue
            diff = float(abs(value - benchmark))
            assertmsg = "Test day, variable %s, %s: %f vs %f" % \
                        (day, var_to_benchmark, benchmark, value)
            self.assertAlmostEqual(diff, 0., precision, assertmsg)

    def runTest(self):
        for eachvar, precision in self.benchmark_vars:
            # Skip testing rooting depth in potential model because of
            # fortran WOFOST does not limit rooting depth in potential mode.
            if (eachvar == "rd") and (self.mode == "pp"):
                continue
            self.run_benchmark(eachvar, precision)

class TestWaterlimitedSpringBarley1965(PyWofostTestingTemplate):
    crop = 3
    year = 1965     # CAMPAIGNYEAR
    grid = 1        # iso 35042
    mode = "wlp"
    multi_layer = True

class TestPotentialWinterWheat(PyWofostTestingTemplate):
    crop = 1
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedWinterWheat(PyWofostTestingTemplate):
    crop = 1
    year = 2000
    grid = 31031
    mode = "wlp"
        
class TestPotentialGrainMaize(PyWofostTestingTemplate):
    crop = 2
    year = 2000
    grid = 31031
    mode = "pp"
        
class TestWaterlimitedGrainMaize(PyWofostTestingTemplate):
    crop = 2
    year = 2000
    grid = 31031
    mode = "wlp"

class TestPotentialSpringBarley(PyWofostTestingTemplate):
    crop = 3
    year = 2000
    grid = 31031
    mode = "pp"

class TestWaterlimitedSpringBarley(PyWofostTestingTemplate):
    crop = 3
    year = 2000
    grid = 31031
    mode = "wlp"

class TestPotentialPotato(PyWofostTestingTemplate):
    crop = 7
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedPotato(PyWofostTestingTemplate):
    crop = 7
    year = 2000
    grid = 31031
    mode = "wlp"

class TestPotentialWinterRapeseed(PyWofostTestingTemplate):
    crop = 10
    year = 2000
    grid = 31031
    mode = "pp"
        
class TestWaterlimitedWinterRapeseed(PyWofostTestingTemplate):
    crop = 10
    year = 2000
    grid = 31031
    mode = "wlp"

class TestPotentialSunflower(PyWofostTestingTemplate):
    crop = 11
    year = 2000
    grid = 31031
    mode = "pp"
        

class TestWaterlimitedSunflower(PyWofostTestingTemplate):
    crop = 11
    year = 2000
    grid = 31031
    mode = "wlp"

def suite(dsn=None):
    """returns the unit tests for the PyWOFOST model.
    
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
    #tests = [TestPotentialWinterWheat('runTest', dsn),
    #         TestWaterlimitedWinterWheat('runTest', dsn),
    #         TestPotentialGrainMaize('runTest', dsn),
    #         TestPotentialSpringBarley('runTest', dsn),
    #         TestPotentialPotato('runTest', dsn),
    #         TestPotentialWinterRapeseed('runTest', dsn),
    #         TestPotentialSunflower('runTest', dsn)]
    #tests = [TestWaterlimitedSpringBarley1965('runTest', dsn)]
    #tests = [TestWaterlimitedWinterWheat('runTest', dsn)]

    # Shuffle tests in random order in order to test whether different results
    # are obtained when tests are run in a different order.
    random.shuffle(tests)
    suite.addTests(tests)
    return suite
