# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Base classes for creating PCSE simulation units.

In general these classes are not to be used directly, but are to be subclassed
when creating PCSE simulation units.
"""
import types
import logging

from ..traitlets import (HasTraits, List, Float, Int, Instance, Dict, Bool, All)
from .dispatcher import DispatcherObject
from .simulationobject import SimulationObject


class BaseEngine(HasTraits, DispatcherObject):
    """Base Class for Engine to inherit from
    """

    def __init__(self):
        HasTraits.__init__(self)
        DispatcherObject.__init__(self)

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)

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

    @property
    def subSimObjects(self):
        """ Find SimulationObjects embedded within self.
        """

        subSimObjects = []
        defined_traits = self.__dict__["_trait_values"]
        for attr in defined_traits.values():
            if isinstance(attr, SimulationObject):
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

    def zerofy(self):
        """Zerofy the value of all rate variables of any sub-SimulationObjects.
        """
        # Walk over possible sub-simulation objects.
        if self.subSimObjects is not None:
            for simobj in self.subSimObjects:
                simobj.zerofy()
