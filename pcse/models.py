import sqlalchemy as sa

from .engine import Engine
from .traitlets import Int, Enum, Unicode
from .util import merge_dict
from . import exceptions as exc


class _Wofost71Base(Engine):
    """Base class for running WOFOST7.1 simulations

    This class inherits from engine and provides the Engine with the
    configuration file for WOFOST 7.1 Potential Production or water-limited
    production depending on the setting of `self.simulation_mode`. The
    latter is defined below in the subclasses `WOFOST71_PP` and
    `WOFOST71_WLP_FD`

    Moreover, it provides the methods `store_to_database` for sending
    the WOFOST simulation results to the table 'sim_results_timeseries'
    """
    # Definition of run identifiers
    crop_name = Unicode()
    year = Int(-99)
    member_id = 0
    simulation_mode = Enum(["wlp", "pp"])

    def __init__(self, sitedata, timerdata, soildata, cropdata,
                 weatherdataprovider):
        Engine.__init__(self, sitedata, timerdata, soildata, cropdata,
                        weatherdataprovider, config=self.config)

        # Run descriptions
        self.crop_name = cropdata["CRPNAM"]
        self.year = timerdata["CAMPAIGNYEAR"]

    #---------------------------------------------------------------------------
    def store_to_database(self, metadata=None, runid=None):
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
            msg = ("Keyword 'runid' should provide the database columns "+
                   "describing the WOFOST run.")
            raise exc.PCSEError(msg)
        # Merge records with model state variables with the PCSE run ID
        recs = [merge_dict(rec, runid) for rec in self._saved_variables]

        if isinstance(metadata, sa.schema.MetaData):
            table_sim_results_ts = sa.Table('sim_results_timeseries', metadata,
                                            autoload=True)
            i = table_sim_results_ts.insert()
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
        tmsg = "PCSE/WOFOST71_PP OUTPUT FILE FOR: \n" +\
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
        if len(self._saved_variables) == 0:
            print "No simulation results available yet."
            return

        tmsg = "TIME-SERIES RESULTS:\n"
        header = "%10s" % "day"
        for varname in self.mconf.OUTPUT_VARS:
            header += (",%10s" % varname)
        header += "\n"
        tmsg += header
        msg += tmsg

        for rec in self._saved_variables:
            tmsg = "%8s" % rec["day"]
            for varname in self.mconf.OUTPUT_VARS:
                value = rec[varname]
                if value is not None:
                    tmsg += (",%10.3f" % value)
                else:
                    tmsg += (",%10s" % "N/A")
            tmsg += "\n"
            msg += tmsg

        with open(outputfile, 'w') as fp:
            fp.write(msg)


class Wofost71_PP(_Wofost71Base):
    """Convenience class for running WOFOST7.1.

    This class inherits from `_WOFOST71Base` and only sets the flag
    for the potential production level: 'pp' and the right
    configuration file.
    """
    simulation_mode = "pp"
    config = "Wofost71_PP.conf"

class Wofost71_WLP_FD(_Wofost71Base):
    """Convenience class for running WOFOST7.1.

    This class inherits from `_WOFOST71Base` and only sets the flag
    for the production level 'pp' or 'wlp' for water-limited production
    for free draining soils and the right configuration file.
    """
    simulation_mode = "wlp"
    config = "Wofost71_WLP_FD.conf"
