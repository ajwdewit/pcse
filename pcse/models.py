# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from .engine import Engine


class Wofost71_PP(Engine):
    """Convenience class for running WOFOST7.1 Potential Production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost71_PP.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost71_WLP_FD(Engine):
    """Convenience class for running WOFOST7.1 water-limited production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost71_WLP_FD.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINTUL3(Engine):
    """The LINTUL model (Light INTerception and UtiLisation) is a simple general crop model,
    which simulates dry matter production as the result of light interception and utilization
    with a constant light use efficiency.

    LINTUL3 simulates crop growth under water-limited and nitrogen-limited conditions

    :param parameterprovider: A `ParameterProvider` object providing model
        parameters as key/value pairs. The parameterprovider encapsulates
        the different parameter sets for crop, soil and site parameters.
    :param weatherdataprovider: An instance of a WeatherDataProvider that can
        return weather data in a WeatherDataContainer for a given date.
    :param agromanagement: AgroManagement data. The data format is described
        in the section on agronomic management.
    """
    config = "Lintul3.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)