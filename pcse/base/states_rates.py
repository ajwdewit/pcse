# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import logging
from datetime import date

from ..traitlets import (HasTraits, List, Float, Int, Instance, Dict, Bool, All)
from ..pydispatch import dispatcher
from ..util import Afgen
from .. import exceptions as exc
from ..settings import settings
from .variablekiosk import VariableKiosk


class ParamTemplate(HasTraits):
    """Template for storing parameter values.

    This is meant to be subclassed by the actual class where the parameters
    are defined.

    example::

        >>> import pcse
        >>> from pcse.base import ParamTemplate
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
          File "pcse/base.py", line 205, in __init__
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

        valid = lambda s: not (s.startswith("_") or s.startswith("trait"))
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
                self._kiosk.register_variable(id(self), attr, type=self._vartype,
                                              publish=True)
                self.observe(handler=self._update_kiosk, names=attr, type=All)
            else:
                self._kiosk.register_variable(id(self), attr, type=self._vartype,
                                              publish=False)
        # Check if the set of published variables is exhausted, otherwise
        # raise an error.
        if len(publish) > 0:
            msg = ("Unknown variable(s) specified with the publish " +
                   "keyword: %s") % publish
            raise exc.PCSEError(msg)

    # def __setattr__(self, attr, value):
    #     # Attributes starting with "_" can be assigned or updated regardless
    #     # of whether the object is locked.
    #     #
    #     # Note that the check on startswith("_") *MUST* be the first otherwise
    #     # the assignment of some trait internals will fail
    #     if attr.startswith("_"):
    #         HasTraits.__setattr__(self, attr, value)
    #     elif attr in self._valid_vars:
    #         if not self._locked:
    #             HasTraits.__setattr__(self, attr, value)
    #         else:
    #             msg = "Assignment to locked attribute '%s' prevented." % attr
    #             raise AttributeError(msg)
    #     else:
    #         msg = "Assignment to non-existing attribute '%s' prevented." % attr
    #         raise AttributeError(msg)

    def _update_kiosk(self, change):
        """Update the variable_kiosk through trait notification.
        """
        self._kiosk.set_variable(id(self), change["name"], change["new"])

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

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)


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
        >>> from pcse.base import VariableKiosk, StatesTemplate
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
          File "pcse/base.py", line 396, in __init__
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
            msg = ("Initial value given for unknown state variable(s): " +
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
        if name in self.rates:
            # known attribute: set value:
            self.rates[name] = value
        elif not self.__initialized:
            # new attribute: allow whe not yet initialized:
            object.__setattr__(self, name, value)
        else:
            # new attribute: disallow according ancestorial ruls:
            super(StatesWithImplicitRatesTemplate, self).__setattr__(name, value)

    def __getattr__(self, name):
        if name in self.rates:
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
        zero_value = {Bool: False, Int: 0, Float: 0.}

        d = {}
        for name, value in self.traits().items():
            if name not in self._valid_vars:
                continue
            try:
                d[name] = zero_value[value.__class__]
            except KeyError:
                msg = ("Rate variable '%s' not of type Float, Bool or Int. " +
                       "Its zero value cannot be determined and it will " +
                       "not be treated by zerofy().") % name
                self.logger.warning(msg)
        return d

    def zerofy(self):
        """Sets the values of all rate values to zero (Int, Float)
        or False (Boolean).
        """
        self._trait_values.update(self._rate_vars_zero)
