# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
from ..pydispatch import dispatcher


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
