# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl) and Herman Berghuijs (herman.berghuijs@wur.nl), January 2024

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


class Wofost72_WLP_CWB(Engine):
    """Convenience class for running WOFOST7.2 water-limited production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost72_WLP_CWB.conf"

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
Wofost71_WLP_FD = Wofost72_WLP_FD = Wofost72_WLP_CWB


class Wofost73_PP(Engine):
    """Convenience class for running WOFOST7.3 Potential Production.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost73_PP.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost73_WLP_CWB(Engine):
    """Convenience class for running WOFOST7.3 Water=limited Production using the Classic Waterbalance.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost73_WLP_CWB.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost73_WLP_MLWB(Engine):
    """Convenience class for running WOFOST7.3 Water=limited Production using the Multi-layer Waterbalance.

    :param parameterprovider: A ParameterProvider instance providing all parameter values
    :param weatherdataprovider: A WeatherDataProvider object
    :param agromanagement: Agromanagement data
    """
    config = "Wofost73_WLP_MLWB.conf"

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
    config = "Lingra_PP.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINGRA_WLP_FD(Engine):
    config = "Lingra_WLP_FD.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class LINGRA_NWLP_FD(Engine):
    config = "Lingra_NWLP_FD.conf"
    __version__ = "1.0.0"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost81_PP(Engine):
    """Convenience class for running WOFOST8.1 potential production
    """
    config = "Wofost81_PP.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost81_WLP_CWB(Engine):
    """Convenience class for running WOFOST8.1  water-limited production using the classic
    waterbalance.
    """
    config = "Wofost81_WLP_CWB.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)


class Wofost81_WLP_MLWB(Engine):
    """Convenience class for running WOFOST8.1  water-limited production using the multi-layer
    waterbalance.
    """
    config = "Wofost81_WLP_MLWB.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost81_NWLP_CWB_CNB(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the classic waterbalance and classic nitrogen balance.
    """
    config = "Wofost81_NWLP_CWB_CNB.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost81_NWLP_MLWB_CNB(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the multi-layer waterbalance and classic nitrogen balance.
    """
    config = "Wofost81_NWLP_MLWB_CNB.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)

class Wofost81_NWLP_MLWB_SNOMIN(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the multi-layer waterbalance and the SNOMIN carbon/nitrogen balance.
    """
    config = "Wofost81_NWLP_MLWB_SNOMIN.conf"

    def __init__(self, parameterprovider, weatherdataprovider, agromanagement):
        Engine.__init__(self, parameterprovider, weatherdataprovider, agromanagement,
                        config=self.config)