# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Base classes for creating PCSE simulation units.

In general these classes are not to be used directly, but are to be subclassed
when creating PCSE simulation units.
"""
import types
import logging
from datetime import date
import cPickle
from collections import Counter

from .traitlets import (HasTraits, Any, Float, Int, Instance, Dict, Bool,
                        Enum, AfgenTrait)
from .pydispatch import dispatcher
from .util import Afgen
from . import exceptions as exc
from .decorators import prepare_states
from .settings import settings

class VariableKiosk(dict):
    """VariableKiosk for registering and publishing state variables in PCSE.
    
    No parameters are needed for instantiating the VariableKiosk.
    All variables that are
    defined within PCSE will be registered within the VariableKiosk, while
    usually only a small subset of those will be published with the kiosk.
    The value of the published 
    variables can be retrieved with the bracket notation as the variableKiosk
    is essentially a (somewhat fancy) dictionary. 
    
    Registering/deregistering rate and state variables goes through the
    `self.register_variable()` and `self.deregister_variable()` methods while the
    `set_variable()` method is used to update a value of a published variable.
    In general, none of these methods need to be called by users directly as
    the logic within the `StatesTemplate` and `RatesTemplate` takes care of
    this.
    
    Finally, the `variable_exists()` can be used to check if a variable is
    registered, while the `flush_states()` and `flush_rates()` are used to
    remove (flush) the values of any published state and rate variables.
    
    example::

        >>> import pcse
        >>> from pcse.base_classes import VariableKiosk
        >>> 
        >>> v = VariableKiosk()
        >>> id0 = 0
        >>> v.register_variable(id0, "VAR1", type="S", publish=True)
        >>> v.register_variable(id0, "VAR2", type="S", publish=False)
        >>> 
        >>> id1 = 1
        >>> v.register_variable(id1, "VAR3", type="R", publish=True)
        >>> v.register_variable(id1, "VAR4", type="R", publish=False)
        >>> 
        >>> v.set_variable(id0, "VAR1", 1.35)
        >>> v.set_variable(id1, "VAR3", 310.56)
        >>> 
        >>> print v
        Contents of VariableKiosk:
         * Registered state variables: 2
         * Published state variables: 1 with values:
          - variable VAR1, value: 1.35
         * Registered rate variables: 2
         * Published rate variables: 1 with values:
          - variable VAR3, value: 310.56

        >>> print v["VAR3"]
        310.56
        >>> v.set_variable(id0, "VAR3", 750.12)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "pcse/base_classes.py", line 148, in set_variable
            raise exc.VariableKioskError(msg % varname)
        pcse.exceptions.VariableKioskError: Unregistered object tried to set the value of variable 'VAR3': access denied.
        >>> 
        >>> v.flush_rates()
        >>> print v
        Contents of VariableKiosk:
         * Registered state variables: 2
         * Published state variables: 1 with values:
          - variable VAR1, value: 1.35
         * Registered rate variables: 2
         * Published rate variables: 1 with values:
          - variable VAR3, value: undefined
        
        >>> v.flush_states()
        >>> print v
        Contents of VariableKiosk:
         * Registered state variables: 2
         * Published state variables: 1 with values:
          - variable VAR1, value: undefined
         * Registered rate variables: 2
         * Published rate variables: 1 with values:
          - variable VAR3, value: undefined
    """
    
    def __init__(self):
        dict.__init__(self)
        self.registered_states = {}
        self.registered_rates  = {}
        self.published_states = {}
        self.published_rates  = {}
    
    def __setitem__(self, item, value):
        msg = "See set_variable() for setting a variable."
        raise RuntimeError(msg)
    
    def __contains__(self, item):
        """Checks if item is in self.registered_states or self.registered_rates.
        """
        return dict.__contains__(self, item)

    def __str__(self):
        msg = "Contents of VariableKiosk:\n"
        msg += " * Registered state variables: %i\n" % len(self.registered_states)
        msg += " * Published state variables: %i with values:\n" % len(self.published_states)
        for varname in self.published_states:
            if varname in self:
                value = self[varname]
            else:
                value = "undefined"
            msg += "  - variable %s, value: %s\n" % (varname, value)
        msg += " * Registered rate variables: %i\n" % len(self.registered_rates)
        msg += " * Published rate variables: %i with values:\n" % len(self.published_rates)
        for varname in self.published_rates:
            if varname in self:
                value = self[varname]
            else:
                value = "undefined"
            msg += "  - variable %s, value: %s\n" % (varname, value)
        return msg
        
    def register_variable(self, oid, varname, type, publish=False):
        """Register a varname from object with id, with given type
        
        :param oid: Object id (from python builtin id() function) of the
            state/rate object registering this variable.
        :param varname: Name of the variable to be registered, e.g. "DVS" 
        :param type: Either "R" (rate) or "S" (state) variable, is handled
            automatically by the states/rates template class.
        :param publish: True if variable should be published in the kiosk,
            defaults to False
        """

        self._check_duplicate_variable(varname)
        if type.upper() == "R":           
            self.registered_rates[varname] = oid
            if publish is True:
                self.published_rates[varname] = oid
        elif type.upper() == "S":
            self.registered_states[varname] = oid
            if publish is True:
                self.published_states[varname] = oid
        else:
            msg = "Variable type should be 'S'|'R'"
            raise exc.VariableKioskError(msg)

    def deregister_variable(self, oid, varname):
        """Object with id(object) asks to deregister varname from kiosk
        
        :param oid: Object id (from python builtin id() function) of the
            state/rate object registering this variable.
        :param varname: Name of the variable to be registered, e.g. "DVS" 
        """
        if varname in self.registered_states:
            #print "Deregistering '%s'" % varname
            if oid != self.registered_states[varname]:
                msg = "Wrong object tried to deregister variable '%s'." \
                      % varname
                raise exc.VariableKioskError(msg)
            else:
                self.registered_states.pop(varname)
            if varname in self.published_states:
                self.published_states.pop(varname)
        elif varname in self.registered_rates:
            #print "Deregistering '%s'" % varname
            if oid != self.registered_rates[varname]:
                msg = "Wrong object tried to deregister variable '%s'." \
                      % varname
                raise exc.VariableKioskError(msg)
            else:
                self.registered_rates.pop(varname)
            if varname in self.published_rates:
                self.published_rates.pop(varname)
        else:
            msg = "Failed to deregister variabe '%s'!" % varname
            raise exc.VariableKioskError(msg)

        # Finally remove the value from the internal dictionary
        if varname in self:
            self.pop(varname)

    def _check_duplicate_variable(self, varname):
        """Checks if variables are not registered twice.
        """
        if varname in self.registered_rates or \
           varname in self.registered_states:
            msg = "Duplicate state/rate variable '%s' encountered!"
            raise exc.VariableKioskError(msg % varname)
        
    def set_variable(self, id, varname, value):
        """Let object with id, set the value of variable varname

        :param id: Object id (from python builtin id() function) of the
            state/rate object registering this variable.
        :param varname: Name of the variable to be updated
        :param value: Value to be assigned to the variable.       
        """
        
        if varname in self.published_rates:
            if self.published_rates[varname] == id:
                dict.__setitem__(self, varname, value)
            else:
                msg = "Unregistered object tried to set the value "+\
                      "of variable '%s': access denied."
                raise exc.VariableKioskError(msg % varname)
        elif varname in self.published_states:
            if self.published_states[varname] == id:
                dict.__setitem__(self, varname, value)
            else:
                msg = "Unregistered object tried to set the value of variable "+\
                "%s: access denied."
                raise exc.VariableKioskError(msg % varname)
        else:
            msg = "Variable '%s' not published in VariableKiosk."
            raise exc.VariableKioskError(msg % varname)
    
    def variable_exists(self, varname):
        """ Returns True if the state/rate variable is registered in the kiosk.

        :param varname: Name of the variable to be checked for registration.
        """
        
        if varname in self.registered_rates or \
           varname in self.registered_states:
            return True
        else:
            return False
        
    def flush_rates(self):
        """flush the values of all published rate variable from the kiosk.
        """
        for key in self.published_rates.keys():
            self.pop(key, None)

    def flush_states(self):
        """flush the values of all state variable from the kiosk.
        """
        for key in self.published_states.keys():
            self.pop(key, None)


class ParamTemplate(HasTraits):
    """Template for storing parameter values.
    
    This is meant to be subclassed by the actual class where the parameters
    are defined.
    
    example::

        >>> import pcse
        >>> from pcse.base_classes import ParamTemplate
        >>> from pcse.traitlets import Float
        >>> 
        >>> 
        >>> class Parameters(ParamTemplate):
        ...     A = Float()
        ...     B = Float()
        ...     C = Float()
        ... 
        >>> parvalues = {"A" :1., "B" :-99, "C":2.45}
        >>> params = Parameters(parvalues)
        >>> params.A
        1.0
        >>> params.A; params.B; params.C
        1.0
        -99.0
        2.4500000000000002
        >>> parvalues = {"A" :1., "B" :-99}
        >>> params = Parameters(parvalues)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "pcse/base_classes.py", line 205, in __init__
            raise exc.ParameterError(msg)
        pcse.exceptions.ParameterError: Value for parameter C missing.
    """

    def __init__(self, parvalues):
        
        HasTraits.__init__(self)
        
        for parname in self.trait_names():
            # If the attribute of the class starts with "trait" than
            # this is a special attribute and not a WOFOST parameter
            if parname.startswith("trait"):
                continue
            # else check if the parname is available in the dictionary
            # of parvalues
            if parname not in parvalues:
                msg = "Value for parameter %s missing." % parname
                raise exc.ParameterError(msg)
            value = parvalues[parname]
            #if isinstance(getattr(self, parname), (Float, Int, Bool, Enum)):
            #    # Single value parameter
            #    setattr(self, parname, value)
            #else:
            #    # AFGEN table parameter
            #    setattr(self, parname, afgen(value))
            if isinstance(getattr(self, parname), (Afgen)):
                # AFGEN table parameter
                setattr(self, parname, Afgen(value))
            else:
                # Single value parameter
                setattr(self, parname, value)
                
    def __setattr__(self, attr, value):
        if attr.startswith("_"):
            HasTraits.__setattr__(self, attr, value)            
        elif hasattr(self, attr):
            HasTraits.__setattr__(self, attr, value)
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)


#-------------------------------------------------------------------------------
def check_publish(publish):
    """ Convert the list of published variables to a set with unique elements.
    """

    if publish is None:
        publish = []
    elif isinstance(publish, str):
        publish = [publish]      
    elif isinstance(publish, (list, tuple)):
        pass
    else:
        msg = "The publish keyword should specify a string or a list of strings"
        raise RuntimeError(msg)
    return set(publish)
    
#-------------------------------------------------------------------------------
class StatesRatesCommon(HasTraits):
    _kiosk = Instance(VariableKiosk)
    _valid_vars = Instance(set)
    _locked = Bool(False)

    def __init__(self, kiosk=None, publish=None):
        """Set up the common stuff for the states and rates template
        including variables that have to be published in the kiosk
        """
        
        HasTraits.__init__(self)

        # Make sure that the variable kiosk is provided
        if not isinstance(kiosk, VariableKiosk):
            msg = ("Variable Kiosk must be provided when instantiating rate " +
                   "or state variables.")
            raise RuntimeError(msg)
        self._kiosk = kiosk
        
        # Check publish variable for correct usage
        publish = check_publish(publish)
        
        # Determine the rate/state attributes defined by the user
        self._valid_vars = self._find_valid_variables()
        
        # Register all variables with the kiosk and optionally publish them.
        self._register_with_kiosk(publish)

    def _find_valid_variables(self):
        """Returns a set with the valid state/rate variables names. Valid rate
        variables have names not starting with 'trait' or '_'. 
        """
        
        valid = lambda s : not (s.startswith("_") or s.startswith("trait"))
        r = [name for name in self.trait_names() if valid(name)]
        return set(r)
    
    def _register_with_kiosk(self, publish):
        """Register the variable with the variable kiosk.

        Here several operations are carried out:
         1. Register the  variable with the kiosk, if rates/states are
            registered twice an error will be raised, this ensures
            uniqueness of rate/state variables across the entire model.
         2 If the  variable name is included in the list set by publish
           keyword then set a trigger on that variable to update its value 
           in the kiosk.
        
         Note that self._vartype determines if the variables is registered
         as a state variable (_vartype=="S") or rate variable (_vartype=="R")
        """

        for attr in self._valid_vars:
            if attr in publish:
                publish.remove(attr)
                self._kiosk.register_variable(id(self), attr,
                                              type=self._vartype,
                                              publish=True)
                self.on_trait_change(self._update_kiosk, attr)
            else:
                self._kiosk.register_variable(id(self), attr,
                                              type=self._vartype,
                                              publish=False)                    
        
        # Check if the set of published variables is exhausted, otherwise
        # raise an error.
        if len(publish) > 0:
            msg = ("Unknown variable(s) specified with the publish " +
                   "keyword: %s") % publish
            raise exc.PCSEError(msg)

    def __setattr__(self, attr, value):
        # Attributes starting with "_" can be assigned or updated regardless
        # of whether the object is locked.
        #
        # Note that the check on startswith("_") *MUST* be the first otherwise
        # the assignment of some trait internals will fail
        if attr.startswith("_"):
            HasTraits.__setattr__(self, attr, value)
        elif attr in self._valid_vars:
            if not self._locked:
                HasTraits.__setattr__(self, attr, value)
            else:
                msg = "Assignment to locked attribute '%s' prevented." % attr
                raise AttributeError(msg)
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)

    def _update_kiosk(self, trait_name, oldvalue, newvalue):
        """Update the variable_kiosk through trait notification.
        """
        #print "Updating published variable '%s' from %s to %s" % \
        #      (trait_name, oldvalue, newvalue)
        self._kiosk.set_variable(id(self), trait_name, newvalue)
    
    def unlock(self):
        "Unlocks the attributes of this class."
        self._locked = False

    def lock(self):
        "Locks the attributes of this class."
        self._locked = True

    def _delete(self):
        """Deregister the variables from the kiosk before garbage
        collecting.
        
        This method is coded as _delete() and must by explicitly called
        because of precarious handling of __del__() in python.
        """
        for attr in self._valid_vars:
            self._kiosk.deregister_variable(id(self), attr)


#-------------------------------------------------------------------------------
class StatesTemplate(StatesRatesCommon):
    """Takes care of assigning initial values to state variables, registering
    variables in the kiosk and monitoring assignments to variables that are
    published.
        
    :param kiosk: Instance of the VariableKiosk class. All state variables
        will be registered in the kiosk in order to enfore that variable names
        are unique across the model. Moreover, the value of variables that
        are published will be available through the VariableKiosk. 
    :param publish: Lists the variables whose values need to be published 
        in the VariableKiosk. Can be omitted if no variables need to be
        published.

    Initial values for state variables can be specified as keyword when instantiating
    a States class.

    example::

        >>> import pcse
        >>> from pcse.base_classes import VariableKiosk, StatesTemplate
        >>> from pcse.traitlets import Float, Integer, Instance
        >>> from datetime import date
        >>> 
        >>> k = VariableKiosk()
        >>> class StateVariables(StatesTemplate):
        ...     StateA = Float()
        ...     StateB = Integer()
        ...     StateC = Instance(date)
        ... 
        >>> s1 = StateVariables(k, StateA=0., StateB=78, StateC=date(2003,7,3),
        ...                     publish="StateC")
        >>> print s1.StateA, s1.StateB, s1.StateC
        0.0 78 2003-07-03
        >>> print k
        Contents of VariableKiosk:
         * Registered state variables: 3
         * Published state variables: 1 with values:
          - variable StateC, value: 2003-07-03
         * Registered rate variables: 0
         * Published rate variables: 0 with values:
        
        >>> 
        >>> s2 = StateVariables(k, StateA=200., StateB=1240)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "pcse/base_classes.py", line 396, in __init__
            raise exc.PCSEError(msg)
        pcse.exceptions.PCSEError: Initial value for state StateC missing.

    """
    
    _kiosk = Instance(VariableKiosk)
    _locked = Bool(False)
    _vartype = "S"

    def __init__(self, kiosk=None, publish=None, **kwargs):
        
        StatesRatesCommon.__init__(self, kiosk, publish)
        
        # set initial state value
        for attr in self._valid_vars:
            if attr in kwargs:
                value = kwargs.pop(attr)
                setattr(self, attr, value)
            else:
                msg = "Initial value for state %s missing." % attr
                raise exc.PCSEError(msg)
        
        # Check if kwargs is empty, otherwise issue a warning
        if len(kwargs) > 0:
            msg = ("Initial value given for unknown state variable(s): "+
                   "%s") % kwargs.keys()
            logging.warn(msg)
            
        # Lock the object to prevent further changes at this stage.
        self._locked = True
        
        
    def touch(self):
        """Re-assigns the value of each state variable, thereby updating its
        value in the variablekiosk if the variable is published."""
        
        self.unlock()
        for name in self._valid_vars:
            value = getattr(self, name)
            setattr(self, name, value)
        self.lock()
    

#-------------------------------------------------------------------------------
class StatesWithImplicitRatesTemplate(StatesTemplate):
    """Container class for state variables that have an associated rate.

    The rates will be generated upon initialization having the same name as their states,
    prefixed by a lowercase character 'r'.
    After initialization no more attributes can be implicitly added.
    Call integrate() to integrate all states with their current rates; the rates are reset to 0.0.
    
    States are all attributes descending from Float and not prefixed by an underscore.
    """
        
    rates = {}
    __initialized = False

    def __setattr__(self, name, value):
        if self.rates.has_key(name):
            # known attribute: set value:             
            self.rates[name] = value
        elif not self.__initialized:
            # new attribute: allow whe not yet initialized:
            object.__setattr__(self, name, value)
        else:
            # new attribute: disallow according ancestorial ruls:
            super(StatesWithImplicitRatesTemplate, self).__setattr__(name, value)
            
            
    def __getattr__(self, name):
        if self.rates.has_key(name):
            return self.rates[name]
        else:
            object.__getattribute__(self, name)


    def initialize_rates(self):
        self.rates = {}
        self.__initialized = True
        
        for s in self.__class__.listIntegratedStates():
            self.rates['r' + s] = 0.0


    def integrate(self, delta):
        # integrate all:
        for s in self.listIntegratedStates():
            rate = getattr(self, 'r' + s)
            state = getattr(self, s)
            newvalue = state + delta * rate
            setattr(self, s, newvalue)
          
        # reset all rates  
        for r in self.rates:
          self.rates[r] = 0.0
          
          
            
    @classmethod
    def listIntegratedStates(cls):
        return sorted([a for a in cls.__dict__ if isinstance(getattr(cls, a), Float) and not a.startswith('_')])



    @classmethod
    def initialValues(cls):
        return dict((a, 0.0) for a in cls.__dict__ if isinstance(getattr(cls, a), Float) and not a.startswith('_'))





#-------------------------------------------------------------------------------
class RatesTemplate(StatesRatesCommon):
    """Takes care of registering variables in the kiosk and monitoring
    assignments to variables that are published.
        
    :param kiosk: Instance of the VariableKiosk class. All rate variables
        will be registered in the kiosk in order to enfore that variable names
        are unique across the model. Moreover, the value of variables that
        are published will be available through the VariableKiosk. 
    :param publish: Lists the variables whose values need to be published 
        in the VariableKiosk. Can be omitted if no variables need to be
        published.
        
    For an example see the `StatesTemplate`. The only difference is that the
    initial value of rate variables does not need to be specified because
    the value will be set to zero (Int, Float variables) or False (Boolean
    variables).
    """
    
    _rate_vars_zero = Instance(dict)    
    _vartype = "R"

    def __init__(self, kiosk=None, publish=None):
        """Set up the RatesTemplate and set monitoring on variables that
        have to be published.
        """
        
        StatesRatesCommon.__init__(self, kiosk, publish)
        
        # Determine the zero value for all rate variable if possible
        self._rate_vars_zero = self._find_rate_zero_values()
        
        # Initialize all rate variables to zero or False
        self.zerofy()
        
        # Lock the object to prevent further changes at this stage.
        self._locked = True
        
    def _find_rate_zero_values(self):
        """Returns a dict with the names with the valid rate variables names as keys and
        the values are the zero values used by the zerofy() method. This means 0 for Int,
        0.0 for Float en False for Bool.
        """

        # Define the zero value for Float, Int and Bool
        zero_value = {Bool:False, Int:0, Float:0.}
    
        d = {}
        for name, value in self.traits().iteritems():
            if name not in self._valid_vars:
                continue
            try:
                d[name] = zero_value[value.__class__]
            except KeyError:
                msg = ("Rate variable '%s' not of type Float, Bool or Int. "+
                       "Its zero value cannot be determined and it will "+
                       "not be treated by zerofy().") % name
                logging.warn(msg)
        return d
    
    def zerofy(self):
        """Sets the values of all rate values to zero (Int, Float)
        or False (Boolean).
        """
        self._trait_values.update(self._rate_vars_zero)

#-------------------------------------------------------------------------------
class DispatcherObject(object):
    """Class only defines the _send_signal() and _connect_signal() methods.
    
    This class is only to be inherited from, not to be used directly.
    """

    def _send_signal(self, signal, *args, **kwargs):
        """Send <signal> using the dispatcher module.
        
        The VariableKiosk of this SimulationObject is used as the sender of
        the signal. Additional arguments to the _send_signal() method are 
        passed to dispatcher.send()
        """
        
        self.logger.debug("Sent signal: %s" % signal)
        dispatcher.send(signal=signal, sender=self.kiosk, *args, **kwargs)
    
    def _connect_signal(self, handler, signal):
        """Connect the handler to the signal using the dispatcher module.
        
        The handler will only react on signals that have the SimulationObjects
        VariableKiosk as sender. This ensure that different PCSE model instances
        in the same runtime environment will not react to each others signals.
        """
        
        dispatcher.connect(handler, signal, sender=self.kiosk)
        self.logger.debug("Connected handler '%s' to signal '%s'." % (handler, signal))

#-------------------------------------------------------------------------------
class SimulationObject(HasTraits, DispatcherObject):
    """Base class for PCSE simulation objects.
    
    :param day: start date of the simulation
    :param kiosk: variable kiosk of this PCSE instance
    
    The day and kiosk are mandatory variables and must be passed when
    instantiating a SimulationObject. 
    
    """

    # Placeholders for logger, params, states, rates and variable kiosk
    logger = Instance(logging.Logger)
    states = Instance(StatesTemplate)
    rates  = Instance(RatesTemplate)
    params = Instance(ParamTemplate)
    kiosk  = Instance(VariableKiosk)
    
    # Placeholder for a list of sub-SimulationObjects. This is to avoid
    # having to loop through all attributes when doing a variable look-up
    subSimObjects = Instance(list)

    # Placeholder for variables that are to be set during finalizing.
    _for_finalize = Dict()
    
    def __init__(self, day, kiosk, *args, **kwargs):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)

        # Check that day variable is specified
        if not isinstance(day, date):
            msg = ("%s should be instantiated with the simulation start " +
                   "day as first argument!")
            raise exc.PCSEError(msg % loggername)

        # Check that kiosk variable is specified and assign to self
        if not isinstance(kiosk, VariableKiosk):
            msg = ("%s should be instantiated with the VariableKiosk " +
                   "as second argument!")
            raise exc.PCSEError(msg % loggername)
        self.kiosk = kiosk

        self.logger = logging.getLogger(loggername)
        self.initialize(day, kiosk, *args, **kwargs)
        self.subSimObjects = self._find_SubSimObjects()
        self.logger.debug("Component successfully initialized on %s!" % day)

    def initialize(self, *args, **kwargs):
        msg = "`initialize` method not yet implemented on %s" % self.__class__.__name__
        raise NotImplementedError(msg)

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
            rebuild = False
            if isinstance(value, SimulationObject) or \
               isinstance(getattr(self, attr), SimulationObject):
                rebuild = True
            HasTraits.__setattr__(self, attr, value)
            if rebuild is True:
                self.subSimObjects = self._find_SubSimObjects()
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)

    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
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
        if self.subSimObjects is not None:
            while len(self.subSimObjects) > 0:
                obj = self.subSimObjects.pop()
                obj._delete()
            
    #---------------------------------------------------------------------------
    def _find_SubSimObjects(self):
        """ Find SimulationObjects embedded within self.
        """
        
        subSimObjects = []
        defined_traits = self.__dict__["_trait_values"]
        for attr in defined_traits.itervalues():
            if isinstance(attr, SimulationObject):
                #print "Found SimObj: %s" % attr.__class__
                subSimObjects.append(attr)
        return subSimObjects

    #---------------------------------------------------------------------------
    def finalize(self, day):
        """ Run the _finalize call on subsimulation objects
        """
        # Update the states object with the values stored in the _for_finalize
        # dictionary
        if self.states is not None:
            self.states.unlock()
            while len(self._for_finalize) > 0:
                k,v = self._for_finalize.popitem()
                setattr(self.states, k, v)
            self.states.lock()
        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.finalize(day)
                
    #---------------------------------------------------------------------------
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

    #---------------------------------------------------------------------------
    def zerofy(self):
        """Zerofy the value of all rate variables of this and any sub-SimulationObjects.
        """

        if self.rates is not None:
            self.rates.zerofy()

        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.zerofy()

#-------------------------------------------------------------------------------
class AncillaryObject(HasTraits, DispatcherObject):
    """Base class for PCSE ancillary objects.
    
    Ancillary objects do not carry out simulation, but often are useful for
    wrapper objects. Still to have some aspects in common with SimulationObjects
    such as the existence of self.logger and self.kiosk, the locked
    behaviour requiring you to define the class attributes and the possibility
    to send/receive signals.
    """
    
    # Placeholders for logger, variable kiosk and parameters
    logger = Instance(logging.Logger)
    kiosk  = Instance(VariableKiosk)
    params = Instance(ParamTemplate)
    
    #---------------------------------------------------------------------------
    def __init__(self, kiosk, *args, **kwargs):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        self.logger = logging.getLogger(loggername)

        # Check that kiosk variable is specified and assign to self
        if not isinstance(kiosk, VariableKiosk):
            msg = "%s should be instantiated with the VariableKiosk "+\
                  "as second argument!"
            raise RuntimeError(msg % loggername)

        self.kiosk = kiosk
        self.initialize(kiosk, *args, **kwargs)
        self.logger.debug("Component successfully initialized!")

    #---------------------------------------------------------------------------
    def __setattr__(self, attr, value):
        if attr.startswith("_"):
            HasTraits.__setattr__(self, attr, value)            
        elif hasattr(self, attr):
            HasTraits.__setattr__(self, attr, value)
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)


class SlotPickleMixin(object):
    """This mixin makes it possible to pickle/unpickle objects with __slots__ defined.

    In many programs, one or a few classes have a very large number of instances.
    Adding __slots__ to these classes can dramatically reduce the memory footprint
    and improve execution speed by eliminating the instance dictionary. Unfortunately,
    the resulting objects cannot be pickled. This mixin makes such classes pickleable
    again and even maintains compatibility with pickle files created before adding
    __slots__.

    Recipe taken from:
    http://code.activestate.com/recipes/578433-mixin-for-pickling-objects-with-__slots__/
    """
    def __getstate__(self):
        return dict(
            (slot, getattr(self, slot))
            for slot in self.__slots__
            if hasattr(self, slot)
        )

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)


class WeatherDataContainer(SlotPickleMixin):
    """Class for storing weather data elements.

    Weather data elements are provided through keywords that are also the
    attribute names under which the variables can accessed in the
    WeatherDataContainer. So the keyword TMAX=15 sets an attribute
    TMAX with value 15.

    The following keywords are compulsory:

    :keyword LAT: Latitude of location (decimal degree)
    :keyword LON: Longitude of location (decimal degree)
    :keyword ELEV: Elevation of location (meters)
    :keyword DAY: the day of observation (python datetime.date)
    :keyword IRRAD: Incoming global radiaiton (J/m2/day)
    :keyword TMIN: Daily minimum temperature (Celsius)
    :keyword TMAX: Daily maximum temperature (Celsius)
    :keyword VAP: Daily mean vapour pressure (hPa)
    :keyword RAIN: Daily total rainfall (cm/day)
    :keyword WIND: Daily mean wind speed at 2m height (m/sec)
    :keyword E0: Daily evaporation rate from open water (cm/day)
    :keyword ES0: Daily evaporation rate from bare soil (cm/day)
    :keyword ET0: Daily evapotranspiration rate from reference crop (cm/day)

    There are two optional keywords arguments:

    :keyword TEMP: Daily mean temperature (Celsius), will otherwise be
                   derived from (TMAX+TMIN)/2.
    :keyword SNOWDEPTH: Depth of snow cover (cm)
    """
    sitevar = ["LAT", "LON", "ELEV"]
    required = ["IRRAD", "TMIN", "TMAX", "VAP", "RAIN", "E0", "ES0", "ET0", "WIND"]
    optional = ["SNOWDEPTH", "TEMP", "TMINRA"]
    # In the future __slots__ can be extended or attribute setting can be allowed
    # by add '__dict__' to __slots__.
    __slots__ = sitevar + required + optional + ["DAY"]

    units = {"IRRAD": "J/m2/day", "TMIN": "Celsius", "TMAX": "Celsius", "VAP": "hPa",
             "RAIN": "cm/day", "E0": "cm/day", "ES0": "cm/day", "ET0": "cm/day",
             "LAT": "Degrees", "LON": "Degrees", "ELEV": "m", "SNOWDEPTH": "cm",
             "TEMP": "Celsius", "TMINRA": "Celsius", "WIND": "m/sec"}

    # ranges for meteorological variables
    ranges = {"LAT": (-90., 90.),
              "LON": (-180., 180.),
              "ELEV": (-300, 6000),
              "IRRAD": (0., 40e6),
              "TMIN": (-50., 60.),
              "TMAX": (-50., 60.),
              "VAP": (0.06, 199.3),  # hPa, computed as sat. vapour pressure at -50, 60 Celsius
              "RAIN": (0, 25),
              "E0": (0., 2.5),
              "ES0": (0., 2.5),
              "ET0": (0., 2.5),
              "WIND": (0., 100.),
              "SNOWDEPTH": (0., 250.),
              "TEMP": (-50., 60.),
              "TMINRA": (-50., 60.)}

    def __init__(self, *args, **kwargs):

        # only keyword parameters should be used for weather data container
        if len(args) > 0:
            msg = ("WeatherDataContainer should be initialized by providing weather " +
                   "variables through keywords only. Got '%s' instead.")
            raise exc.PCSEError(msg % args)

        # First assign site variables
        for varname in self.sitevar:
            try:
                setattr(self, varname, float(kwargs.pop(varname)))
            except (KeyError, ValueError) as e:
                msg = "Site parameter '%s' missing or invalid when building WeatherDataContainer: %s"
                raise exc.PCSEError(msg, varname, e)

        # check if we have a DAY element
        if "DAY" not in kwargs:
            msg = "Date of observations 'DAY' not provided when building WeatherDataContainer."
            raise exc.PCSEError(msg)
        self.DAY = kwargs.pop("DAY")

        # Loop over required arguments to see if all required variables are there
        for varname in self.required:
            value = kwargs.pop(varname, None)
            try:
                setattr(self, varname, float(value))
            except (KeyError, ValueError, TypeError) as e:
                msg = "%s: Weather attribute '%s' missing or invalid numerical value: %s"
                logging.warning(msg, self.DAY, varname, value)

        # Loop over optional arguments
        for varname in self.optional:
            value = kwargs.pop(varname, None)
            if value is None:
                continue
            else:
                try:
                    setattr(self, varname, float(value))
                except (KeyError, ValueError, TypeError) as e:
                    msg = "%s: Weather attribute '%s' missing or invalid numerical value: %s"
                    logging.warning(msg, self.DAY, varname, value)

        # Check for remaining unknown arguments
        if len(kwargs) > 0:
            msg = "WeatherDataContainer: unknown keywords '%s' are ignored!"
            logging.warning(msg, kwargs.keys())

    def __setattr__(self, key, value):
        # Override to allow range checking on known meteo variables.
        if key in self.ranges:
            vmin, vmax = self.ranges[key]
            if not vmin <= value <= vmax:
                msg = "Value (%s) for meteo variable '%s' outside allowed range (%s, %s)." % (value, key, vmin, vmax)
                raise exc.PCSEError(msg)
        SlotPickleMixin.__setattr__(self, key, value)

    def __str__(self):
        msg = "Weather data for %s (DAY)\n" % self.DAY
        for v in self.required:
            value = getattr(self, v, None)
            if value is None:
                msg += "%5s: element missing!\n"
            else:
                unit = self.units[v]
                msg += "%5s: %12.2f %9s\n" % (v, value, unit)
        for v in self.optional:
            value = getattr(self, v, None)
            if value is None:
                continue
            else:
                unit = self.units[v]
                msg += "%5s: %12.2f %9s\n" % (v, value, unit)
        msg += ("Latitude  (LAT): %8.2f degr.\n" % self.LAT)
        msg += ("Longitude (LON): %8.2f degr.\n" % self.LON)
        msg += ("Elevation (ELEV): %6.1f m.\n" % self.ELEV)
        return msg
                
    def add_variable(self, varname, value, unit):
        """Adds an attribute <varname> with <value> and given <unit>
        
        :param varname: Name of variable to be set as attribute name (string)
        :param value: value of variable (attribute) to be added.
        :param unit: string representation of the unit of the variable. Is
            only use for printing the contents of the WeatherDataContainer.
        """
        if varname not in self.units:
            self.units[varname] = unit
        setattr(self, varname, value)

#-------------------------------------------------------------------------------
class WeatherDataProvider(object):
    """Base class for all weather data providers.
    
    Support for weather ensembles in a WeatherDataProvider has to be indicated
    by setting the class variable `supports_ensembles = True`
    
    Example::
    
        class MyWeatherDataProviderWithEnsembles(WeatherDataProvider):
            supports_ensembles = True
            
            def __init__(self):
                WeatherDataProvider.__init__(self)

                # remaining initialization stuff goes here.
    """
    supports_ensembles = False

    # Descriptive items for a WeatherDataProvider
    longitude = None
    latitude = None
    elevation = None
    description = None
    _first_date = None
    _last_date = None
    angstA = None
    angstB = None
    # model used for reference ET
    ETmodel = "PM"

    def __init__(self):
        self.store = {}
        # Define a logger
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        self.logger = logging.getLogger(loggername)

    def _dump(self, cache_fname):
        """Dumps the contents into cache_fname using cPickle.

        Dumps the values of self.store, longitude, latitude, elevation and description
        """
        with open(cache_fname, "wb") as fp:
            dmp = (self.store, self.elevation, self.longitude, self.latitude, self.description, self.ETmodel)
            cPickle.dump(dmp, fp, cPickle.HIGHEST_PROTOCOL)

    def _load(self, cache_fname):
        """Loads the contents from cache_fname using cPickle.

        Loads the values of self.store, longitude, latitude, elevation and description
        from cache_fname and also sets the self.first_date, self.last_date
        """

        with open(cache_fname, "rb") as fp:
            (store, self.elevation, self.longitude, self.latitude, self.description, ETModel) = cPickle.load(fp)

        # Check if the reference ET from the cache file is calculated with the same model as
        # specified by self.ETmodel
        if ETModel != self.ETmodel:
            msg = "Mismatch in reference ET from cache file."
            raise exc.PCSEError(msg)

        self.store.update(store)

    @property
    def first_date(self):
        try:
            self._first_date = min(self.store)[0]
        except ValueError:
            pass
        return self._first_date

    @property
    def last_date(self):
        try:
            self._last_date = max(self.store)[0]
        except ValueError:
            pass
        return self._last_date

    @property
    def missing(self):
        missing = (self.last_date - self.first_date).days - len(self.store) + 1
        return missing

    def check_keydate(self, key):
        """Check representations of date for storage/retrieval of weather data.
        
        The following formats are supported:
        
        1. a date object
        2. a datetime object 
        3. a string of the format YYYYMMDD
        4. a string of the format YYYYDDD
        
        Formats 2-4 are all converted into a date object internally.
        """

        import datetime as dt
        if isinstance(key, dt.datetime):
            return key.date()
        elif isinstance(key, dt.date):
            return key
        elif isinstance(key, str):
            skey = key.strip()
            l = len(skey)
            if l==8:
                # assume YYYYMMDD
                dkey = dt.datetime.strptime(skey,"%Y%m%d")
                return dkey.date()
            elif l==7:
                # assume YYYYDDD
                dkey = dt.datetime.strptime(skey,"%Y%j")
                return dkey.date()
            else:
                msg = "Key for WeatherDataProvider not recognized as date: %s"
                raise KeyError(msg % key)
        else:
            msg = "Key for WeatherDataProvider not recognized as date: %s"
            raise KeyError(msg % key)
    
    def _store_WeatherDataContainer(self, wdc, keydate, member_id=0):
        """Stores the WDC under given keydate and member_id.
        """
        
        if member_id != 0 and self.supports_ensembles is False:
            msg = "Storing ensemble weather is not supported."
            raise exc.WeatherDataProviderError(msg)

        kd = self.check_keydate(keydate)
        if not (isinstance(member_id, int) and member_id >= 0):
            msg = "Member id should be a positive integer, found %s" % member_id
            raise exc.WeatherDataProviderError(msg)

        self.store[(kd, member_id)] = wdc
    
    def __call__(self, day, member_id=0):
        
        if self.supports_ensembles is False and member_id != 0:
            msg = "Retrieving ensemble weather is not supported by %s" % self.__class__.__name__
            raise exc.WeatherDataProviderError(msg)

        keydate = self.check_keydate(day)
        
        if self.supports_ensembles is False:
            if self.logger is not None:
                msg = "Retrieving weather data for day %s" % keydate
                self.logger.debug(msg)
            try:
                return self.store[(keydate, 0)]
            except KeyError as e:
                msg = "No weather data for %s." % keydate
                raise exc.WeatherDataProviderError(msg)
        else:
            if self.logger is not None:
                msg = "Retrieving ensemble weather data for day %s member %i" % \
                      (keydate, member_id)
                self.logger.debug(msg)
            try:
                return self.store[(keydate, member_id)]
            except KeyError:
                msg = "No weather data for (%s, %i)." % (keydate,member_id)
                raise exc.WeatherDataProviderError(msg)

    def __str__(self):

        msg = "Weather data provided by: %s\n" % self.__class__.__name__
        msg += "--------Description---------\n"
        if isinstance(self.description, str):
            msg += ("%s\n" % self.description)
        else:
            for l in self.description:
                msg += ("%s\n" % str(l))
        msg += "----Site characteristics----\n"
        msg += "Elevation: %6.1f\n" % self.elevation
        msg += "Latitude:  %6.3f\n" % self.latitude
        msg += "Longitude: %6.3f\n" % self.longitude
        msg += "Data available for %s - %s\n" % (self.first_date, self.last_date)
        msg += "Number of missing days: %i\n" % self.missing
        return msg

class BaseEngine(HasTraits, DispatcherObject):
    """Base Class for Engine to inherit from
    """
    # Placeholders for logger, params, states, rates and variable kiosk
    logger = Instance(logging.Logger)

    # Placeholder for a list of sub-SimulationObjects. This is to avoid
    # having to loop through all attributes when doing a variable look-up
    subSimObjects = Instance(list)

    def __init__(self):
        HasTraits.__init__(self)
        DispatcherObject.__init__(self)

        # Define logger
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        self.logger = logging.getLogger(loggername)

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
            rebuild = False
            if isinstance(value, SimulationObject) or \
               isinstance(getattr(self, attr), SimulationObject):
                rebuild = True
            HasTraits.__setattr__(self, attr, value)
            if rebuild is True:
                self.subSimObjects = self._find_SubSimObjects()
        else:
            msg = "Assignment to non-existing attribute '%s' prevented." % attr
            raise AttributeError(msg)

    def _find_SubSimObjects(self):
        """ Find SimulationObjects embedded within self.
        """

        subSimObjects = []
        defined_traits = self.__dict__["_trait_values"]
        for attr in defined_traits.itervalues():
            if isinstance(attr, SimulationObject):
                #print "Found SimObj: %s" % attr.__class__
                subSimObjects.append(attr)
        return subSimObjects

    def get_variable(self, varname):
        """ Return the value of the specified state or rate variable.

        :param varname: Name of the variable.

        Note that the `get_variable()` will first search for `varname` exactly
        as specified (case sensitive). If the variable cannot be found, it will
        look for the uppercase name of that variable. This is purely for
        convenience.
        """

        # Check if variable is registered in the kiosk, also check for
        # name in upper case as most variables are defined in upper case.
        # If variable is not registered in the kiosk then return None directly.
        if self.kiosk.variable_exists(varname):
            v = varname
        elif self.kiosk.variable_exists(varname.upper()):
            v = varname.upper()
        else:
            return None

        if v in self.kiosk:
            return self.kiosk[v]

        # Search for variable by traversing the hierarchy
        value = None
        for simobj in self.subSimObjects:
            value = simobj.get_variable(v)
            if value is not None:
                break
        return value

    #---------------------------------------------------------------------------
    def zerofy(self):
        """Zerofy the value of all rate variables of any sub-SimulationObjects.
        """

        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.zerofy()


class ParameterProvider(HasTraits):
    """Simple class providing a dictionary-like single interface for parameter values.

    The idea behind this class is twofold. First of all by encapsulating the four
    different parameter types (e.g. sitedata, timerdata, etc) into a single object,
    the signature of the `initialize()` method of each `SimulationObject` can be
    harmonized across all SimulationObjects. Second, the ParameterProvider itself
    can be easily adapted when different sets of parameter values are needed. For
    example when running PCSE with crop rotations, different sets of timerdata and
    cropdata are needed, this can now be handled easily by enhancing
    ParameterProvider to rotate new sets of timerdata and cropdata on a CROP_FINISH
    signal.

    See also the `MultiCropParameterProvider`
    """
    _maps = list()
    _sitedata = dict()
    _soildata = dict()
    _cropdata = dict()
    _timerdata = dict()
    _override = dict()

    def __init__(self, sitedata=None, timerdata=None, soildata=None, cropdata=None):
        if sitedata is not None:
            self._sitedata = sitedata
        else:
            self._sitedata = {}
        if cropdata is not None:
            self._cropdata = cropdata
        else:
            self._cropdata = {}
        if soildata is not None:
            self._soildata = soildata
        else:
            self._soildata = {}
        if timerdata is not None:
            self._timerdata = timerdata
        else:
            self._timerdata = {}
        self._override = {}
        self._maps = [self._override, self._sitedata, self._timerdata, self._soildata, self._cropdata]
        self._test_uniqueness()

    def set_crop_type(self, crop_id=None, crop_start_type=None, crop_end_type=None):
        """Set the start_type and end type of the crop which is relevant for
        the phenology module.

        :param crop_id: string identifying the crop type, is ignored as only
               one crop is assumed to be here.
        :param crop_start_type: start type for the given crop: 'sowing'|'emergence'
        :param crop_end_type: end type for the given crop: 'maturity'|'harvest'|'earliest'
        """

        self._timerdata["CROP_START_TYPE"] = crop_start_type
        self._timerdata["CROP_END_TYPE"] = crop_end_type
        self._test_uniqueness()

    def set_override(self, varname, value, check=True):
        """"Override the value of variable varname in the parameterprovider.

        Overriding the value of particular variable is often useful for example
        when running for different sets of parameters or for calibration
        purposes.

        Note that if check=True (default) varname should already exist in one of site, timer,
        soil or cropdata.
        """

        if check:
            if varname in self:
                self._override[varname] = value
            else:
                msg = "Cannot override '%s', parameter does not exist." % varname
                raise exc.PCSEError(msg)
        else:
            self._override[varname] = value

    def clear_override(self, varname=None):
        """Removes variable varname from override, without arguments all variables are removed."""

        if varname is None:
            self._override.clear()
        else:
            if varname in self._override:
                self._override.pop(varname)
            else:
                msg = "Cannot clear varname '%s' from override" % varname

    def _test_uniqueness(self):
        # Check if parameter names are unique
        parnames = []
        for mapping in [self._sitedata, self._timerdata, self._soildata, self._cropdata]:
            parnames.extend(mapping.keys())
        unique = Counter(parnames)
        for parname, count in unique.items():
            if count > 1:
                msg = "Duplicate parameter found: %s" % parname
                raise RuntimeError(msg)

    def __getitem__(self, key):
        for mapping in self._maps:
            try:
                return mapping[key]
            except KeyError:
                pass
        raise KeyError(key)

    def __contains__(self, key):
        for mapping in self._maps:
            if key in mapping:
                return True
        return False

class MultiCropParameterProvider(ParameterProvider):
    """Parameter provider that allows multiple crop
    parameter sets to be specified. This ParameterProvider is
    designed to be combined with the AgroManager in order to
    facilitate crop rotations.

    Note that timerdata does not have to be provided anymore
    because this role has been taken over by the AgroManager.

    :param sitedata: A dictionary with site parameters
    :param soildata: A dictionary with soil parameters
    :param multi_cropdata: A dict of dicts with the crop parameters
        keyed on the `crop_id` used in the crop calendar. For
        example:  multi_cropdata = {'winter-wheat': {<parameters for winter-wheat>},
                                    'maize': {<parameters for maize>},
                                    etc...
    :keyword max_root_depth_name: The name of the crop parameter specifying the
        maximum crop rooting depth. Defaults to "RDMCR"
    :keyword init_root_depth_name: The name of the crop parameter specifying the
        initial crop rooting depth. Defaults to "RDI".

    warning:: The WOFOST ClassicWaterBalance needs to know the maximum rooting depth
    of all crops (RDMCR) in order to define its soil layers. Therefore all
    sets of crop parameters are searched for the maximum value of RDMCR. Moreover, the
    initial rooting depth (RDI) needs to be the same across all crop types. The names
    of these parameter can be provided through the keywords specified above but default
    to 'RDMCR' and 'RDI'.
    """
    _multi_cropdata = dict()
    _RDMCR_max = 0.  # maximum rooting depth across all crop types for the water balance
    _RDI = 0.   # Initial rooting depth

    def __init__(self, sitedata, soildata, multi_cropdata, max_root_depth_name="RDMCR",
                 init_root_depth_name="RDI"):

        self._sitedata = sitedata
        self._soildata = soildata
        self._cropdata = {}
        self._timerdata = {}
        self._multi_cropdata = multi_cropdata

        # Get maximum rooting depth and initial rooting depth over all crops
        RDI = []
        for crop_id, cropdata in self._multi_cropdata.items():
            if cropdata[max_root_depth_name] > self._RDMCR_max:
                self._RDMCR_max = cropdata[max_root_depth_name]
            RDI.append(cropdata[init_root_depth_name])

        # Test if all crops have the same initial rooting depth
        if len(set(RDI)) > 1:
            msg = "Initial rooting depth (%s) not the same across all crop types." % init_root_depth_name
            raise exc.PCSEError(msg)
        self._RDI = RDI[0]

        # update the cropdata to provide default values for RDI and RDMCR which
        # are need by the waterbalance at the initialization.
        self._cropdata[max_root_depth_name] = self._RDMCR_max
        self._cropdata[init_root_depth_name] = self._RDI

        self._maps = [self._sitedata, self._timerdata, self._soildata, self._cropdata]
        self._test_uniqueness()

    def set_crop_type(self, crop_id=None, crop_start_type=None, crop_end_type=None):
        """Switch the crop parameters in the MultiCropParameterProvider to crop type given by crop_id.

        :param crop_id: string identifying the crop type
        :param crop_start_type: start type for the given crop: 'sowing'|'emergence'
        :param crop_end_type: end type for the given crop: 'maturity'|'harvest'|'earliest'
        """

        if crop_id not in self._multi_cropdata:
            msg = "Crop parameters for crop (%s) cannot be found in the multi_cropdata." % crop_id
            raise exc.PCSEError(msg)

        self._cropdata.clear()
        self._cropdata.update(self._multi_cropdata[crop_id])
        self._timerdata["CROP_START_TYPE"] = crop_start_type
        self._timerdata["CROP_END_TYPE"] = crop_end_type

        self._maps = [self._sitedata, self._timerdata, self._soildata, self._cropdata]
        self._test_uniqueness()
