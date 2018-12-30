# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Base classes for creating PCSE simulation units.

In general these classes are not to be used directly, but are to be subclassed
when creating PCSE simulation units.
"""
import types
import logging
from datetime import date
import pickle
from collections import Counter, MutableMapping

from ..traitlets import (HasTraits, List, Float, Int, Instance, Dict, Bool, All)
from ..pydispatch import dispatcher
from ..util import Afgen
from .. import exceptions as exc
from ..settings import settings
from .variablekiosk import VariableKiosk


class ParameterProvider(MutableMapping):
    """Class providing a dictionary-like interface over all parameter sets (crop, soil, site).
    It acts very much like a ChainMap with some additional features.

    The idea behind this class is threefold. First of all by encapsulating the
    different parameter sets (sitedata, cropdata, soildata) into a single object,
    the signature of the `initialize()` method of each `SimulationObject` can be
    harmonized across all SimulationObjects. Second, the ParameterProvider itself
    can be easily adapted when different sets of parameter values are needed. For
    example when running PCSE with crop rotations, different sets of cropdata
    are needed, this can now be handled easily by enhancing
    ParameterProvider to rotate a new set of cropdata when the engine receives a
    CROP_START signal. Finally, specific parameter values can be easily changed
    by setting an `override` on that parameter.

    See also the `MultiCropDataProvider`
    """
    _maps = list()
    _sitedata = dict()
    _soildata = dict()
    _cropdata = dict()
    _timerdata = dict()
    _override = dict()
    _unique_parameters = list()
    _iter = 0  # Counter for iterator
    _ncrops_activated = 0  # Counts the number of times `set_crop_type()` has been called.

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

    def set_active_crop(self, crop_name=None, variety_name=None, crop_start_type=None, crop_end_type=None):
        """Activate the crop parameters for the given crop_name and variety_name.

        :param crop_name: string identifying the crop name, is ignored as only
               one crop is assumed to be here.
        :param variety_name: string identifying the variety name, is ignored as only
               one crop is assumed to be here.
        :param crop_start_type: start type for the given crop: 'sowing'|'emergence'
        :param crop_end_type: end type for the given crop: 'maturity'|'harvest'|'earliest'

        In case of crop rotations, there is a new set of crop parameters needed when a new
        crop is started. This routine activates the crop parameters for the given crop_name and
        variety_name. The `crop_name`, `variety_name` `crop_start_type` and `crop_end_type`
        are defined in the agromanagement and supported by the AgroManager.

        Note that many CropDataProviders are not designed for crop rotations and only support a single
        crop whose parameters are active by default. In this case a call to `set_active_crop()` has no
        effect and the `crop_name` and `variety_name` parameters are ignored.
        CropDataProviders that support crop rotations explicitly have to subclass from
        `pcse.base.MultiCropDataProvider` in order to be recognized.

        Besides the crop parameters, this method also sets the `crop_start_type` and `crop_end_type` of the
        crop which is required for all crops by the phenology module.

        """

        self._timerdata["CROP_START_TYPE"] = crop_start_type
        self._timerdata["CROP_END_TYPE"] = crop_end_type
        if isinstance(self._cropdata, MultiCropDataProvider):
            # we have a MultiCropDataProvider, so set the active crop and variety
            self._cropdata.set_active_crop(crop_name, variety_name)
        else:
            # we do not have a MultiCropDataProvider, this means that crop rotations are not supported
            # At the first call this is OK. However issue a warning with subsequent calls
            # to set_crop_type() are done because we cannot change the set of crop parameters
            if self._ncrops_activated == 0:
                pass
            else:
                # has been called multiple times
                msg = "A second crop was scheduled: however, the CropDataProvider does not " \
                      "support multiple crop parameter sets. This will only work for crop" \
                      "rotations with the same crop."
                self.logger.warning(msg)

        self._ncrops_activated += 1
        self._test_uniqueness()

    @property
    def logger(self):
        loggername = "%s.%s" % (self.__class__.__module__,
                                self.__class__.__name__)
        return logging.getLogger(loggername)

    def set_override(self, varname, value, check=True):
        """"Override the value of parameter varname in the parameterprovider.

        Overriding the value of particular parameter is often useful for example
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
        """Removes parameter varname from the set of overridden parameters.

        Without arguments all overridden parameters are removed.
        """

        if varname is None:
            self._override.clear()
        else:
            if varname in self._override:
                self._override.pop(varname)
            else:
                msg = "Cannot clear varname '%s' from override" % varname
                raise exc.PCSEError(msg)

    def _test_uniqueness(self):
        """Check if parameter names are unique and raise an error if duplicates occur.

        Note that the uniqueness is not tested for parameters in self._override as this
        is specifically meant for overriding parameters.
        """
        parnames = []
        for mapping in [self._sitedata, self._timerdata, self._soildata, self._cropdata]:
            parnames.extend(mapping.keys())
        unique = Counter(parnames)
        for parname, count in unique.items():
            if count > 1:
                msg = "Duplicate parameter found: %s" % parname
                raise exc.PCSEError(msg)

    @property
    def _unique_parameters(self):
        """Returns a list of unique parameter names across all sets of parameters.

        This includes the parameters in self._override in order to be able to
        iterate over all parameters in the ParameterProvider.
        """
        s = []
        for mapping in self._maps:
            s.extend(mapping.keys())
        return sorted(list(set(s)))

    def __getitem__(self, key):
        """Returns the value of the given parameter (key).

        Note that the search order in self._map is such that self._override is tested first for the
        existence of the key. Thus ensuring that overridden parameters will be found first.

        :param key: parameter name to return
        """
        for mapping in self._maps:
            if key in mapping:
                return mapping[key]
        raise KeyError(key)

    def __contains__(self, key):
        for mapping in self._maps:
            if key in mapping:
                return True
        return False

    def __str__(self):
        msg = "ParameterProvider providing %i parameters, %i parameters overridden: %s."
        return msg % (len(self), len(self._override), self._override.keys())

    def __setitem__(self, key, value):
        """Override an existing parameter (key) by value.

         The parameter that is overridden is added to self._override, note that only *existing*
         parameters may be overridden this way. If it is needed to really add a *new* parameter
         than use: ParameterProvider.set_override(key, value, check=False)

        :param key: The name of the parameter to override
        :param value: the value of the parameter
        """
        if key in self:
            self._override[key] = value
        else:
            msg = "Cannot override parameter '%s', parameter does not exist. " \
                  "to bypass this check use: set_override(parameter, value, check=False)" % key
            raise exc.PCSEError(msg)

    def __delitem__(self, key):
        """Deletes a parameter from self._override.

        Note that only parameters that exist in self._override can be deleted. This also means that
        if an parameter is overridden its original value will return after a parameter is deleted.

        :param key: The name of the parameter to delete
        """
        if key in self._override:
            self._override.pop(key)
        elif key in self:
            msg = "Cannot delete default parameter: %s" % key
            raise exc.PCSEError(msg)
        else:
            msg = "Parameter not found!"
            raise KeyError(msg)

    def __len__(self):
        return len(self._unique_parameters)

    def __iter__(self):
        return self

    def next(self):
        i = self._iter
        if i < len(self):
            self._iter += 1
            return self._unique_parameters[self._iter - 1]
        else:
            self._iter = 0
            raise StopIteration


class MultiCropDataProvider(dict):

    def __init__(self):
        dict.__init__(self)
        self._store = {}

    def set_active_crop(self, crop_name, variety_name):
        """Sets the crop parameters for the crop identified by crop_name and variety_name.

        Needs to be implemented by each subclass of MultiCropDataProvider
        """
        msg = "'set_crop_type' method should be implemented specifically for each" \
              "subclass of MultiCropDataProvider."
        raise NotImplementedError(msg)

