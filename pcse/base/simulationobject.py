# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import types
import logging
from datetime import date

from .dispatcher import DispatcherObject
from ..traitlets import (HasTraits, List, Float, Int, Instance, Dict, Bool, All)
from .. import exceptions as exc
from .variablekiosk import VariableKiosk
from .states_rates import StatesTemplate, RatesTemplate, ParamTemplate



class SimulationObject(HasTraits, DispatcherObject):
    """Base class for PCSE simulation objects.

    :param day: start date of the simulation
    :param kiosk: variable kiosk of this PCSE instance

    The day and kiosk are mandatory variables and must be passed when
    instantiating a SimulationObject.

    """

    # Placeholders for logger, params, states, rates and variable kiosk
    states = Instance(StatesTemplate)
    rates = Instance(RatesTemplate)
    params = Instance(ParamTemplate)
    kiosk = Instance(VariableKiosk)

    # Placeholder for variables that are to be set during finalizing.
    _for_finalize = Dict()

    def __init__(self, day, kiosk, *args, **kwargs):
        HasTraits.__init__(self, *args, **kwargs)

        # Check that day variable is specified
        if not isinstance(day, date):
            this = "%s.%s" % (self.__class__.__module__, self.__class__.__name__)
            msg = ("%s should be instantiated with the simulation start " +
                   "day as first argument!")
            raise exc.PCSEError(msg % this)

        # Check that kiosk variable is specified and assign to self
        if not isinstance(kiosk, VariableKiosk):
            this = "%s.%s" % (self.__class__.__module__, self.__class__.__name__)
            msg = ("%s should be instantiated with the VariableKiosk " +
                   "as second argument!")
            raise exc.PCSEError(msg % this)
        self.kiosk = kiosk

        self.initialize(day, kiosk, *args, **kwargs)
        self.logger.debug("Component successfully initialized on %s!" % day)

    def initialize(self, *args, **kwargs):
        msg = "`initialize` method not yet implemented on %s" % self.__class__.__name__
        raise NotImplementedError(msg)

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)

    def integrate(self, *args, **kwargs):
        msg = "`integrate` method not yet implemented on %s" % self.__class__.__name__
        raise NotImplementedError(msg)

    def calc_rates(self, *args, **kwargs):
        msg = "`calc_rates` method not yet implemented on %s" % self.__class__.__name__
        raise NotImplementedError(msg)

    def __setattr__(self, attr, value):
        # __setattr__ has been modified  to enforce that class attributes
        # must be defined before they can be assigned. There are a few
        # exceptions:
        # 1 if an attribute name starts with '_'  it will be assigned directly.
        # 2 if the attribute value is a  function (e.g. types.FunctionType) it
        #   will be assigned directly. This is needed because the
        #   'prepare_states' and 'prepare_rates' decorators assign the wrapped
        #   functions 'calc_rates', 'integrate' and optionally 'finalize' to
        #   the Simulation Object. This will collide with __setattr__ because
        #   these class methods are not defined attributes.
        #
        # Finally, if the value assigned to an attribute is a SimulationObject
        #   or if the existing attribute value is a SimulationObject than
        #   rebuild the list of sub-SimulationObjects.

        if attr.startswith("_") or type(value) is types.FunctionType:
            HasTraits.__setattr__(self, attr, value)
        elif hasattr(self, attr):
            HasTraits.__setattr__(self, attr, value)
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)

    def get_variable(self, varname):
        """ Return the value of the specified state or rate variable.

        :param varname: Name of the variable.

        Note that the `get_variable()` will searches for `varname` exactly
        as specified (case sensitive).
        """

        # Search for variable in the current object, then traverse the hierarchy
        value = None
        if hasattr(self.states, varname):
            value = getattr(self.states, varname)
        elif hasattr(self.rates, varname):
            value = getattr(self.rates, varname)
        # Query individual sub-SimObject for existence of variable v
        else:
            for simobj in self.subSimObjects:
                value = simobj.get_variable(varname)
                if value is not None:
                    break
        return value

    def set_variable(self, varname, value, incr):
        """ Sets the value of the specified state or rate variable.

        :param varname: Name of the variable to be updated (string).
        :param value: Value that it should be updated to (float)
        :param incr: dict that will receive the increments to the updated state
            variables.

        :returns: either the increment of the variable (new - old) or `None`
          if the call was unsuccessful in finding the class method (see below).

        Note that 'setting'  a variable (e.g. updating a model state) is much more
        complex than just `getting` a variable, because often some other
        internal variables (checksums, related state variables) must be updated
        as well. As there is no generic rule to 'set' a variable it is up to
        the model designer to implement the appropriate code to do the update.

        The implementation of `set_variable()` works as follows. First it will
        recursively search for a class method on the simulationobjects with the
        name `_set_variable_<varname>` (case sensitive). If the method is found,
        it will be called by providing the value as input.

        So for updating the crop leaf area index (varname 'LAI') to value '5.0',
        the call will be: `set_variable('LAI', 5.0)`. Internally, this call will
        search for a class method `_set_variable_LAI` which will be executed
        with the value '5.0' as input.
        """
        method_name = "_set_variable_%s" % varname.strip()
        try:
            method_obj = getattr(self, method_name)
            rv = method_obj(value)
            if not isinstance(rv, dict):
                msg = ("Method %s on '%s' should return a dict with the increment of the " +
                       "updated state variables!") % (method_name, self.__class__.__name__)
                raise exc.PCSEError(msg)
            incr.update(rv)
        except AttributeError:  # method is not present: just continue
            pass
        except TypeError:  # method is present but is not callable: error!
            msg = ("Method '%s' on '%s' could not be called by 'set_variable()': " +
                   "check your code!") % (method_name, self.__class__.__name__)
            raise exc.PCSEError(msg)

        for simobj in self.subSimObjects:
            simobj.set_variable(varname, value, incr)

    def _delete(self):
        """ Runs the _delete() methods on the states/rates objects and recurses
        trough the list of subSimObjects.
        """
        if self.states is not None:
            self.states._delete()
            self.states = None
        if self.rates is not None:
            self.rates._delete()
            self.rates = None
        for obj in self.subSimObjects:
            obj._delete()

    @property
    def subSimObjects(self):
        """ Return SimulationObjects embedded within self.
        """

        subSimObjects = []
        defined_traits = self.__dict__["_trait_values"]
        for attr in defined_traits.values():
            if isinstance(attr, SimulationObject):
                subSimObjects.append(attr)
        return subSimObjects

    def finalize(self, day):
        """ Run the _finalize call on subsimulation objects
        """
        # Update the states object with the values stored in the _for_finalize dictionary
        if self.states is not None:
            self.states.unlock()
            while len(self._for_finalize) > 0:
                k, v = self._for_finalize.popitem()
                setattr(self.states, k, v)
            self.states.lock()
        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.finalize(day)

    def touch(self):
        """'Touch' all state variables of this and any sub-SimulationObjects.

        The name comes from the UNIX `touch` command which does nothing on the
        contents of a file but only updates the file metadata (time, etc).
        Similarly, the `touch` method re-assigns the state of each state
        variable causing any triggers (e.g. `on_trait_change()`) to go off.
        This will guarantee that these state values remain available in the
        VariableKiosk.
        """

        if self.states is not None:
            self.states.touch()
        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.touch()

    def zerofy(self):
        """Zerofy the value of all rate variables of this and any sub-SimulationObjects.
        """

        if self.rates is not None:
            self.rates.zerofy()

        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.zerofy()


class AncillaryObject(HasTraits, DispatcherObject):
    """Base class for PCSE ancillary objects.

    Ancillary objects do not carry out simulation, but often are useful for
    wrapper objects. Still to have some aspects in common with SimulationObjects
    such as the existence of self.logger and self.kiosk, the locked
    behaviour requiring you to define the class attributes and the possibility
    to send/receive signals.
    """

    # Placeholders for logger, variable kiosk and parameters
    kiosk = Instance(VariableKiosk)
    params = Instance(ParamTemplate)

    def __init__(self, kiosk, *args, **kwargs):
        HasTraits.__init__(self, *args, **kwargs)

        # Check that kiosk variable is specified and assign to self
        if not isinstance(kiosk, VariableKiosk):
            this = "%s.%s" % (self.__class__.__module__, self.__class__.__name__)
            msg = "%s should be instantiated with the VariableKiosk " \
                  "as second argument!"
            raise RuntimeError(msg % this)

        self.kiosk = kiosk
        self.initialize(kiosk, *args, **kwargs)
        self.logger.debug("Component successfully initialized!")

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)

    def __setattr__(self, attr, value):
        if attr.startswith("_"):
            HasTraits.__setattr__(self, attr, value)
        elif hasattr(self, attr):
            HasTraits.__setattr__(self, attr, value)
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)
