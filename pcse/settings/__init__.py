# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
import os, sys
import importlib
import inspect

from . import default_settings

class Settings(object):
    """
    Settings for PCSE.

    Default values will be read from the module pcse.settings.default_settings.py
    User settings are read from $HOME/.pcse/user_settings.py and override the default ones;
    see the default settings file for a list of all possible variables.
    """

    def __setattr__(self, name, value):
        if name == "METEO_CACHE_DIR":
            if not os.path.exists(value):
                os.mkdir(value)
        if name == "LOG_DIR":
            if not os.path.exists(value):
                os.mkdir(value)
        object.__setattr__(self, name, value)


    def __init__(self):
        # update this dict from global settings (but only for ALL_CAPS settings)
        for setting in dir(default_settings):
            if setting.isupper():
                setattr(self, setting, getattr(default_settings, setting))
            elif setting.startswith("_"):
                pass
            else:
                msg = ("Warning: settings should be ALL_CAPS. Setting '%s' in default_" +
                       "settings will be ignored.") % setting
                print(msg)

        try:
            mod = importlib.import_module("user_settings")
        except ImportError as e:
            raise ImportError(
                ("Could not import settings '%s' (Is it on sys.path? Is there an import" +
                 " error in the settings file?): %s") % ("$HOME/.pcse/user_settings.py", e)
            )

        for setting in dir(mod):
            if setting.isupper():
                setattr(self, setting, getattr(mod, setting))
            elif setting.startswith("_"):
                pass
            else:
                msg = ("Warning: settings should be ALL_CAPS. Setting '%s' in user_" +
                       "settings will be ignored.") % setting
                print(msg)

# Initialize the settings from default_settings and users_settings
settings = Settings()
