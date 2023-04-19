# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014

from .engine import Engine


class Wofost72_PP(Engine):
    """Convenience class for running WOFOST7.2 Potential Production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost72_PP.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost72_WLP_FD(Engine):
    """Convenience class for running WOFOST7.2 water-limited production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost72_WLP_FD.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost72_Phenology(Engine):
    """Convenience class for running WOFOST7.2 phenology only.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost72_Pheno.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


# This is to ensure that old code keeps working
Wofost71_PP = Wofost72_PP
Wofost71_WLP_FD = Wofost72_WLP_FD


class Wofost81_PP(Engine):
    """Convenience class for running WOFOST8.1 potential production (includes N dynamics)

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost81_PP.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost81_WLP_FD(Engine):
    """Convenience class for running WOFOST8.0 water-limited production (includes NPK dynamics)

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost81_WLP_FD.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost81_NWLP_FD(Engine):
    """Convenience class for running WOFOST8.1 nitrogen and water-limited production

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost81_NWLP_FD.conf"

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


class FAO_WRSI(Engine):
    """Convenience class for computing actual crop water use using the Water Requirements
    Satisfaction Index with a (modified) FAO WRSI approach.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """

    config = "FAO_WRSI.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINGRA_PP(Engine):
    """An implementation of the LINGRA model for potential production conditions.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Lingra_PP.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINGRA_WLP_FD(Engine):
    """An implementation of the LINGRA model for water-limited production conditions.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Lingra_WLP_FD.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINGRA_NWLP_FD(Engine):
    """An implementation of the LINGRA model for water and nitrogen limited production conditions.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Lingra_NWLP_FD.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class ALCEPAS(Engine):
    """ALCEPAS Onion growth model for potential production conditions.

    For the model description, see:

    ALCEPAS, an onion growth model based on SUCROS87.I. Development of the model.
     C.L.M. De Visser. 1993. https://doi.org/10.1080/14620316.1994.11516482
    """
    config = "alcepas.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

