# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from .. import exceptions as exc


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
        >>> from pcse.base import VariableKiosk
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
          File "pcse/base.py", line 148, in set_variable
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
        self.registered_rates = {}
        self.published_states = {}
        self.published_rates = {}

    def __setitem__(self, item, value):
        msg = "See set_variable() for setting a variable."
        raise RuntimeError(msg)

    def __contains__(self, item):
        """Checks if item is in self.registered_states or self.registered_rates.
        """
        return dict.__contains__(self, item)

    def __getattr__(self, item):
        """Allow use of attribute notation (eg "kiosk.LAI") on published rates or states.
        """
        return dict.__getitem__(self, item)

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
            # print "Deregistering '%s'" % varname
            if oid != self.registered_states[varname]:
                msg = "Wrong object tried to deregister variable '%s'." \
                      % varname
                raise exc.VariableKioskError(msg)
            else:
                self.registered_states.pop(varname)
            if varname in self.published_states:
                self.published_states.pop(varname)
        elif varname in self.registered_rates:
            # print "Deregistering '%s'" % varname
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
                msg = "Unregistered object tried to set the value " + \
                      "of variable '%s': access denied."
                raise exc.VariableKioskError(msg % varname)
        elif varname in self.published_states:
            if self.published_states[varname] == id:
                dict.__setitem__(self, varname, value)
            else:
                msg = "Unregistered object tried to set the value of variable " \
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
