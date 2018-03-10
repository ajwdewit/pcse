# -*- coding: utf-8 -*-
# Copyright (c) 2004-2018 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""
This module is here only to ensure that all PCSE modules can import internally from .traitlets
while this module loads the actual traitlets modules from the correct location. Moreover,
some traits are adapted to allow `None` as default values and coercing value to float().

Currently an adapted version of the traitlets package is used 'traitlets_pcse'. In the future
the default traitlets package may be used when the `observe()` functionality on `type=All` is implemented.
"""
from traitlets_pcse import *
import traitlets_pcse as tr


class Instance(tr.Instance):

    def __init__(self, *args, **kwargs):
        if 'allow_none' not in kwargs:
            kwargs['allow_none'] = True
        tr.Instance.__init__(self, *args, **kwargs)


class Enum(tr.Enum):

    def __init__(self, *args, **kwargs):
        if 'allow_none' not in kwargs:
            kwargs['allow_none'] = True
        tr.Enum.__init__(self, *args, **kwargs)


class Unicode(tr.Unicode):

    def __init__(self, *args, **kwargs):
        if 'allow_none' not in kwargs:
            kwargs['allow_none'] = True
        tr.Unicode.__init__(self, *args, **kwargs)


class Bool(tr.Bool):

    def __init__(self, *args, **kwargs):
        if 'allow_none' not in kwargs:
            kwargs['allow_none'] = True
        tr.Bool.__init__(self, *args, **kwargs)


class Float(tr.Float):

    def __init__(self, *args, **kwargs):
        if 'allow_none' not in kwargs:
            kwargs['allow_none'] = True
        tr.Float.__init__(self, *args, **kwargs)

    def validate(self, obj, value):
        try:
            value = float(value)
        except:
            self.error(obj, value)
        return value
