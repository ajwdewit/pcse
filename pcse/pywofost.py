#!/usr/bin/env python
import datetime
import copy

import sqlalchemy as sa
try:
    from pandas import DataFrame
except ImportError: # not Available, return a lambda doing nothing
    DataFrame = lambda x: x
    
from .traitlets import HasTraits, Float, Int, Instance, Bool, Enum, Unicode
from .base_classes import VariableKiosk, TopPyWOFOSTObject
from .soil.classic_waterbalance import WaterbalancePP, WaterbalanceFD
from .soil.waterbalance import WaterbalanceLayered
from .soilcropsimulation import SoilCropSimulation
from .timer import Timer
from .agromanagement import AgroManagementSingleCrop as AgroManagement
from .util import merge_dict
from . import exceptions as exc

class PyWofost(TopPyWOFOSTObject):
    """PyWOFOST top level object
            
    Note that this top-level object does not carry out any simulations itself.
    It only takes care of the aspects that are specific to PyWOFOST such as
    1) registering the run descriptor (crop_no,grid_no,year,day,simulation_mode,
    member_id) for which the simulation is being carried out and 2) sending the
    simulation results to the corresponding table which uses the run
    descriptor as its primary key.
    
    The actual soil-crop simulation is carried out within the SimulationObject
    `SoilCropSimulation`. The latter is oblivious to details as run descriptors
    and can therefore be easily reused in other top-level objects. Good
    candidates are developing a ClassicWOFOST which reads/writes to files and
    a PyCGMS which is compatible with the CGMS database.
    """
    # Definition of run identifiers
    crop_name = Unicode()
    year      = Int(-99)
    member_id = Int(-99)
    simulation_mode = Enum(["wlp","pp"])
    
    # place holder for the actual model
    soilcropsimulation = Instance(SoilCropSimulation)
            
    def initialize(self, day, sitedata, timerdata, soildata, cropdata,
                    weatherdataprovider, mode='pp', member_id=0, metadata=None):
        """
        :param day: The start date of the simulation
        :param sitedata: A dictionary(-like) object containing key/value pairs with
            parameters that are specific for this site but not related to the crop,
            the crop calendar or the soil. Examples are the initial conditions of
            the waterbalance such as the initial amount of soil moisture and
            surface storage.
        :param timerdata: A dictionary(-like) object containing key/value pairs
            with parameters related to the system start date, crop calendar, start
            type (sowing|emergence) and end type (maturity|harvest|earliest)
        :param soildata: A dictionary(-like) object containing key/value pairs
            with parameters related to soil where the simulation has to be
            performed.
        :param cropdata: A dictionary(-like) object containing key/value pairs with
            WOFOST crop parameters.
        :param weatherdataprovider: An instance of a WeatherDataProvider that can
            return weather data in a WeatherDataContainer
        :keyword mode: The mode of the simulation to be carried out. Either
            potential production 'pp' or water-limited production 'wlp' (free
            drainage). Defaults to 'pp'
        :keyword member_id: The member id of the PyWOFOST instance. This is only
            used for ensemble simulations. Defaults to zero.
        :keyword metadata: an SQLAlchemy metadata object, defaults to `None`. The
            metadata object will be used to inspect the table `pywofost_output`.
            The columns in the table will be mapped to `get_variable()` thus
            additional output variables can be retrieved by simply adding columns
            to the table `pywofost_output`.
        """
        
        # Run descriptions
        self.crop_name = cropdata["CRPNAM"]
        self.year      = timerdata["CAMPAIGNYEAR"]
        self.member_id = member_id
        self.simulation_mode = mode.lower()
        
        msg = ("Starting new simulation for crop %s, year %4i, " +
               "mode '%s'") % (self.crop_name, self.year, self.simulation_mode)
        self.logger.info(msg)
        #print(msg)

        # Variable Kiosk
        self.kiosk = VariableKiosk()
        
        # Setup the timer
        start_date = timerdata["START_DATE"]
        final_date = timerdata["CROP_END_DATE"]
        timer = Timer(start_date, self.kiosk, final_date)

        # Initialize AgroManagement module
        agromanagement = AgroManagement(start_date, self.kiosk, timerdata,
                                        soildata, sitedata, cropdata)

        # Initialize waterbalance
        if mode.lower() == "pp":
            waterbalance = WaterbalancePP(start_date, self.kiosk, soildata)
        elif "NSL" in soildata: # Multi-layer soil water balance
            waterbalance = WaterbalanceLayered(start_date, self.kiosk, \
                cropdata, soildata, sitedata, {"GW":False,"ZTI":20,"DD":-99})
        else: # Single layer soil water balance
            waterbalance = WaterbalanceFD(start_date, self.kiosk, cropdata,
                                              soildata, sitedata)
        
        # Determine which variables to save from the PyWofost database
        varlist = self._register_variables_to_save(metadata)
        
        # initialize crop-soil model
        self.soilcropsimulation = \
            SoilCropSimulation(start_date, self.kiosk, timer, waterbalance,
                               weatherdataprovider, agromanagement,
                               variables_to_save=varlist)
        
    #---------------------------------------------------------------------------
    def store_to_database(self, metadata=None, runid=None):
        """Stores saved variables of the model run in a database table.
        
        :param metadata: An SQLAlchemy metadata object providing access to the
                         database where the table 'pywofost_output' can be
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
            msg = ("Keyword 'runid' should provide the database columns "+
                   "describing the WOFOST run.")
            raise exc.PCSEError(msg)
        # Merge records with model state variables with the PyWOFOST run ID
        recs = [merge_dict(rec, runid) for rec in \
                self.soilcropsimulation._saved_variables]

        if isinstance(metadata, sa.schema.MetaData):
            table_pyw_output = sa.Table('pywofost_output', metadata,
                                        autoload=True)
            i = table_pyw_output.insert()
            i.execute(recs)
        else:
            msg = ("Keyword metadata should provide an SQLAlchemy " +
                   "MetaData object.")
            raise RuntimeError(msg)

    #---------------------------------------------------------------------------
    def store_to_file(self, outputfile=None):
        """Store simulation results to <outputfile>.
        """

        sep = "--------------------------------------------------------------\n"
        tmsg = "PyWOFOST OUTPUT FILE FOR: \n" +\
               "crop name: %s\n" +\
               "year: %4i\n" +\
               "mode: '%s'\n" 
        msg = tmsg % (self.crop_name, self.year, self.simulation_mode)
        msg += sep
        
        # Collect summary results
        tmsg = "SUMMARY RESULTS:\n"
        summary = [("Day of sowing","DOS","%8s"),("Day of emergence","DOE","%8s"),
                   ("Day of anthesis","DOA","%8s"),("Day of maturity","DOM","%8s"),
                   ("Day of harvest","DOH","%8s"),("Maximum LAI","LAIMAX","%5.2f"),
                   ("Total biomass","TAGP","%7.1f"), ("Yield","TWSO","%7.1f"),
                   ("Harvest Index","HI","%7.3f"), ("Total transpiration","CTRAT","%7.3f")]
        for desc, varname, fmt in summary:
            value = self.get_variable(varname)
            if value is not None:
                templ = "%20s: " + fmt + "\n"
                tmsg += (templ % (desc, value))
            else:
                templ = "%20s: N/A\n"
                tmsg += (templ % (desc))
        msg += tmsg
        msg += sep
        
        # Find variables available in daily records
        if len(self.soilcropsimulation._saved_variables)>0:
            firstrec = copy.deepcopy(self.soilcropsimulation._saved_variables[0])
            firstrec.pop("day")
            varnames = sorted(firstrec.keys())
        else:
            print "No simulation results available yet."
            return

        tmsg = "TIME-SERIES RESULTS:\n"
        header = "%10s" % "day"
        for varname in varnames:
            header += (",%10s" % varname)
        header += "\n"
        tmsg += header
        msg += tmsg

        for rec in self.soilcropsimulation._saved_variables:
            tmsg = "%8s" % rec["day"]
            for varname in varnames:
                value = rec[varname]
                if value is not None:
                    tmsg += (",%10.3f" % value)
                else:
                    tmsg += (",%10s" % "N/A")
            tmsg += "\n"
            msg += tmsg
        
        with open(outputfile, 'w') as fp:
            fp.write(msg)

    #---------------------------------------------------------------------------
    def get_variable(self, varname):
        """Returns the value of variable `varname` or `None` if it cannot be
        found.
        """
        value = self.soilcropsimulation.get_variable(varname)
        return value        
    
    #---------------------------------------------------------------------------
    def grow(self, *args, **kwargs):
        "Maps to `run()` only here for compatibility"
        self.run(*args, **kwargs)

    #---------------------------------------------------------------------------
    def run(self, days=1):
        """Advances the model state with given number of days (default days=1)
        """
        self.soilcropsimulation.run(days)

    #---------------------------------------------------------------------------
    def _register_variables_to_save(self, metadata):
        """Returns the variable names of the variables that will be saved
        during model run.
        
        If metadata is None then a default set of variables is returned.
        Otherwise, metadata must be an SQLALchemy metadata object which is
        used to find out which columns are in the 'pywofost_output' table.
        The column names that are not part of the run_descriptors are returned
        as variable names.
        """

        default = ["dvs","lai","tagp", "twso", "twlv", "twst",
                   "twrt", "tra", "rd", "sm", "wwlow"]
        
        if metadata is None:
            return default
        elif isinstance(metadata, sa.schema.MetaData):
            var = []
            pw_output = sa.Table('pywofost_output', metadata, autoload=True)
            for col in pw_output.columns:
                if col.primary_key not in (1,True):
                    var.append(col.name)
            return var
        else:
            msg = ("Metadata passed is not an SQLAlchemy metadata object. " +
                   "Returning default list of variables.")
            warnings.warn(msg)
            self.logger.warn(msg)
            return default
    
    def get_results(self):
        """Returns the time-series of simulation results.
        
        If the `pandas` package is available then PyWOFOST will return the
        results as a pandas DataFrame with an index on the day. If not, you
        will get a list of dictionaries.
        """
        df = DataFrame(self.soilcropsimulation._saved_variables)
        df = df.set_index("day")
        return df