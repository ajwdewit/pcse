#!/usr/bin/env python
from collections import namedtuple, deque
import logging
import datetime

from .pydispatch import dispatcher
from .traitlets import Instance, Bool
from .base_classes import VariableKiosk, SimulationObject, WeatherDataProvider,\
      AncillaryObject, WeatherDataContainer
from .timer import Timer
from . import signals
from . import exceptions as exc

class SoilCropSimulation(SimulationObject):
    """SimulationObject implementing combined soil/crop system.
    
    `SoilCropSimulation` handles the actual simulation of the combined soil-
    crop system. The central part of the  `SoilCropSimulation` is the
    waterbalance which is continuously simulating during the entire run. In
    contrast, `CropSimulation` objects are only initialized after receiving a
    "CROP_START" signal from the AgroManagement unit. From that point onward,
    the combined soil-crop is simulated including the interactions between
    the soil and crop such as root growth and transpiration.
    
    Similarly, the crop simulation is finalized when receiving a "CROP_FINISH"
    signal. At that moment the `finalize()` section on the cropsimulation is
    executed. Moreover, the "CROP_FINISH" signal can specify that the
    cropsimulation object should be deleted from the hierarchy. The latter is
    useful for further extensions of PyWOFOST for running crop rotations.
    
    Finally, the entire simulation is terminated when a "TERMINATE" signal is
    received. At that point, the `finalize()` section on the waterbalance is 
    executed and the simulation stops.

    **Signals handled by SoilCropSimulation:**
    
    `SoilCropSimulation` handles the following signals:
        * CROP_START: Starts an instance of `CropSimulation` for simulating crop
          growth. See the `_on_CROP_START` handler for details.
        * CROP_FINISH: Runs the `finalize()` section an instance of 
          `CropSimulation` and optionally deletes the cropsimulation instance.
          See the `_on_CROP_FINISH` handler for details.
        * TERMINATE: Runs the `finalize()` section on the waterbalance module
          and terminates the entire simulation.
          See the `_on_TERMINATE` handler for details.
        * OUTPUT:  Preserves a copy of the value of selected state/rate 
          variables for later use.
          See the `_on_OUTPUT` handler for details.
    
    """

    # sub components for simulation
    cropsimulation = Instance(SimulationObject)
    waterbalance   = Instance(SimulationObject)
    agromanagement = Instance(AncillaryObject)
    weatherdataprovider   = Instance(WeatherDataProvider)
    timer = Instance(Timer)
    day   = Instance(datetime.date)
    
    # flags that are being set by signals
    flag_terminate    = Bool(False)
    flag_crop_finish  = Bool(False)
    flag_crop_start   = Bool(False)
    flag_crop_delete  = Bool(False)
    flag_output       = Bool(False)
    
    # placeholders for variables saved during model execution
    _saved_variables = Instance(list)
    _variables_to_save = Instance(list)
    
    # Helper variables
    drv = Instance(WeatherDataContainer)
    TMNSAV = Instance(deque)
    
    def initialize(self, day, kiosk, timer, waterbalance, weatherdataprovider,
                    agromanagement, variables_to_save=None):
        """
        :param day: The start date of the simulation
        :param kiosk: `VariableKiosk` instance
        :param timer: `Timer` Instance
        :param waterbalance: Instance of `WaterbalancePP` or `WaterbalanceFD`
        :param weatherdataprovider: Instance of `WeatherDataProvider` for
               retrieving weather data during the simulation.
        :param agromanagement: AgroManagement object for handling AgroManagement
               actions
        :param variables_to_save: list of variable names that should be
               saved when OUTPUT signals are received.
        """

        # register handlers for starting/finishing the crop simulation, for
        # handling output and terminating the system
        self._connect_signal(self._on_CROP_START, signal=signals.crop_start)
        self._connect_signal(self._on_CROP_FINISH, signal=signals.crop_finish)
        self._connect_signal(self._on_OUTPUT, signal=signals.output)
        self._connect_signal(self._on_TERMINATE, signal=signals.terminate)

        # Variable kiosk for registering and publishing variables
        self.kiosk = kiosk

        # Timer and starting day
        self.timer = timer
        self.day = self.timer()

        # Waterbalance
        self.waterbalance = waterbalance

        # Driving variables
        self.weatherdataprovider = weatherdataprovider
        self.drv = self._get_driving_variables(self.day)

        # Call AgroManagement module for management actions at initialization
        self.agromanagement = agromanagement
        self.agromanagement(self.day, self.drv)

        # Register states to be saved during model run and placeholder
        self._variables_to_save = self._register_variables_to_save(variables_to_save)
        self._saved_variables = list()

        #self._do_some_printing('ini')
        
        # Calculate initial rates
        self.calc_rates(self.day, self.drv)
        
    #---------------------------------------------------------------------------
    def calc_rates(self, day, drv):

        # Start rate calculation on individual components
        if self.cropsimulation is not None:
            self.cropsimulation.calc_rates(day, drv)
        if self.waterbalance is not None:
            self.waterbalance.calc_rates(day, drv)

        # Save state variables of the model
        if self.flag_output:
            self._save_variables(day)

        # Check if flag is present to finish crop simulation
        if self.flag_crop_finish:
            self._finish_cropsimulation(day)

    #---------------------------------------------------------------------------
    def integrate(self, day):

        # Flush state variables from the kiosk before state updates
        self.kiosk.flush_states()

        if self.cropsimulation is not None:
            self.cropsimulation.integrate(day)

        self.waterbalance.integrate(day)

        # Set all rate variables to zero
        self.zerofy()

        # Flush rate variables from the kiosk after state updates
        self.kiosk.flush_rates()

    #---------------------------------------------------------------------------
    def grow(self, *args, **kwargs):
        """Calls run() method, just for backward compatibility"""
        self.run(*args, **kwargs)
        
    #---------------------------------------------------------------------------
    def run(self, days=1):
        """Advances the system state with given number of days"""

        days_done = 0
        while (days_done < days) and (self.flag_terminate is False):
            days_done += 1

            # Update timer
            self.day = self.timer()

            # State integration
            self.integrate(self.day)

            # Driving variables
            self.drv = self._get_driving_variables(self.day)

            # Agromanagement decisions
            self.agromanagement(self.day, self.drv)            

            # Rate calculation
            self.calc_rates(self.day, self.drv)
        
        if self.flag_terminate is True:
            self.waterbalance.finalize(self.day)

    #---------------------------------------------------------------------------
    def _on_CROP_FINISH(self, day, crop_delete=False):
        """Sets the variable 'flag_crop_finish' to True when the signal
        CROP_FINISH is received.
        
        The flag is needed because finishing the crop simulation is deferred to
        the correct place in the processing loop and is done by the routine
        _finish_cropsimulation().
        
        If crop_delete=True the CropSimulation object will be deleted from the
        hierarchy in _finish_cropsimulation(). 
        """
        self.flag_crop_finish = True
        self.flag_crop_delete = crop_delete
        
    #---------------------------------------------------------------------------
    def _on_CROP_START(self, day, cropsimulation):
        """Receives a crop simulation object from the agromanagement unit.
        This object is than assigned to self.cropsimulation and will used
        for simulation of crop dynamics.

        """
        self.logger.debug("Received signal 'CROP_START' on day %s" % day)

        if self.cropsimulation is not None:
            msg = "Want to start a new crop while old crop is still present."
            raise exc.PCSEError(msg)
        self.cropsimulation = cropsimulation
    
    #---------------------------------------------------------------------------
    def _on_TERMINATE(self):
        """Sets the variable 'flag_terminate' to True when the signal TERMINATE
        was received.
        """
        self.flag_terminate = True
        
    #---------------------------------------------------------------------------
    def _on_OUTPUT(self):
        """Sets the variable 'flag_output to True' when the signal OUTPUT
        was received.
        """
        self.flag_output = True
        
    #---------------------------------------------------------------------------
    def _finish_cropsimulation(self, day):
        """Finishes the CropSimulation object when variable 'flag_crop_finish'
        has been set to True based on the signal 'CROP_FINISH' being
        received.
        """
        self.flag_crop_finish = False

        # Run the finalize section of the cropsimulation and sub-components
        self.cropsimulation.finalize(day)
        
        #print "Model started  on day %s" % self.get_variable("doe")
        #print "Model Anthesis on day %s" % self.get_variable("doa")
        #print "Model finishes on day %s on %s" % \
        #    (self.get_variable("dof"), self.get_variable("finish"))
        #print "Final biomass: %7.1f" % self.get_variable("tagp")
        #print "Maximum LAI: %6.2f" % self.get_variable("laimax")
        
        #self._do_some_printing('fin')
        
        # Only remove the crop simulation object from the system when the crop
        # is finished, when explicitly asked to do so.
        if self.flag_crop_delete:
            self.flag_crop_delete = False
            self.cropsimulation._delete()
            self.cropsimulation = None

    #---------------------------------------------------------------------------
    def _get_driving_variables(self, day):
        """Get driving variables, compute derived properties and return it.
        """
        drv = self.weatherdataprovider(day)
        
        # average temperature and average daytemperature (Celsius?)
        drv.add_variable("TEMP", (drv.TMIN + drv.TMAX)/2., "Celcius")
        drv.add_variable("DTEMP", (drv.TEMP + drv.TMAX)/2., "Celcius")

        #  7 day running average of minimum temperature
        TMINRA = self._7day_running_avg(drv.TMIN)
        drv.add_variable("TMINRA", TMINRA, "Celsius")
        
        return drv

    #---------------------------------------------------------------------------
    def _7day_running_avg(self, TMIN):
        """Calculate 7-day running mean of minimum temperature.
        """
        # if self.TMNSAV is None, then initialize a deque of size 7
        if self.TMNSAV is None:
            self.TMNSAV = deque(maxlen=7)
        # Append new value
        self.TMNSAV.appendleft(TMIN)

        return sum(self.TMNSAV)/len(self.TMNSAV)
    
    #---------------------------------------------------------------------------
    def _do_some_printing(self, time):
        """temporary for printing FORTRAN style output.
        """
        if time == "ini" or time == "fin":
            print "\n DAY   WLV   WST    WSO   TAGP   LAI  RD    SM RESRV AVAIL RAIN TRA EVA wet dry"
            print "     kg/ha kg/ha  kg/ha  kg/ha m2/m2  cm vol.fr  cm    cm    mm mm/d mm/d  days"
        
        dic = self.kiosk
        print "%s %5.0f. %4.0f. %5.0f. %5.0f. %4.2f %3.0i. %.3f %4.1f %4.1f %3.0f %.2f %5.2f %2i %2i" % \
            (time, dic.get("WLV",-9), dic.get("WST",-9), dic.get("WSO",-9), dic.get("TAGP",-9), \
             dic.get("LAI",-9), dic.get("RD",-9), dic.get("SM",-9), dic.get("WLOW",-9), dic.get("W",-9), \
             10*dic.get("RAINT",-0.9), dic.get("TRA",-9), dic.get("WTRAT",-9),\
             dic.get("IDOST",-9), dic.get("IDWST",-9))
        
        if time == "fin":
            print "\n SUMMARY :                                                        stress days"
            print " HALT ANTH TWRT   TWLV   TWST   TWSO   TAGP HINDEX TRC  GASST  MREST wet dry "
            print " %s %s %5.0f. %5.0f. %5.0f. %5.0f. %5.0f. %5.2f -9.99  %5.0f.  %4.0f.  %2i  %2i" % \
                (time, dic.get("DOA",'--'), dic.get("TWRT",-9), dic.get("TWLV",-9), dic.get("TWST",-9), dic.get("TWSO",-9), \
                 dic.get("TAGP",-9), dic.get("HI",-9.99), dic.get("GASST",-9999), dic.get("MREST",-999), dic.get("IDOST",-9), dic.get("IDWST",-9))
        
    #---------------------------------------------------------------------------
    def _save_variables(self, day):
        """Appends selected model variables to self._saved_variables for this day.
        """
        # Switch off the flag for generating output
        self.flag_output = False
        
        # find current value of variables to are to be saved
        states = {"day":day}
        for var in self._variables_to_save:
            states[var] = self.get_variable(var)
        self._saved_variables.append(states)

    #---------------------------------------------------------------------------
    def _register_variables_to_save(self, variables_to_save):
        """Returns the variable names of the variables that will be saved
        during model run.
        
        If variable_to_save is None then a default set of variables is returned.
        Otherwise, the specified list is returned."""

        default = ["DVS","LAI","TAGP", "TWSO", "TWLV", "TWST",
                   "TWRT", "TRA", "RD", "SM"]
        
        if variables_to_save is None:
            return default
        else:
            return variables_to_save
