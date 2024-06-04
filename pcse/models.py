# -*- coding: utf-8 -*-
# Copyright (c) 2004-2024 Wageningen Environmental Research, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl) and Herman Berghuijs (herman.berghuijs@wur.nl), May 2024

from .engine import Engine


class Wofost72_PP(Engine):
    """Convenience class for running WOFOST7.2 Potential Production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost72_PP.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.2"
    __waterbalance__ = None
    __nitrogenbalance__ = None


class Wofost72_WLP_CWB(Engine):
    """Convenience class for running WOFOST7.2 water-limited production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost72_WLP_CWB.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.2"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = None


class Wofost72_Phenology(Engine):
    """Convenience class for running WOFOST7.2 phenology only.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost72_Pheno.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.2"
    __waterbalance__ = None
    __nitrogenbalance__ = None


class Wofost73_PP(Engine):
    """Convenience class for running WOFOST7.3 Potential Production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost73_PP.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.3"
    __waterbalance__ = None
    __nitrogenbalance__ = None


class Wofost73_WLP_CWB(Engine):
    """Convenience class for running WOFOST7.3 Water=limited Production using the Classic Waterbalance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost73_WLP_CWB.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.3"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = None


class Wofost73_WLP_MLWB(Engine):
    """Convenience class for running WOFOST7.3 Water=limited Production using the Multi-layer Waterbalance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost73_WLP_MLWB.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "7.3"
    __waterbalance__ = "MLWB"
    __nitrogenbalance__ = None


class Lintul10_NWLP_CWB_CNB(Engine):
    """The LINTUL model (Light INTerception and UtiLisation) is a simple general crop model,
    which simulates dry matter production as the result of light interception and utilization
    with a constant light use efficiency.

    In literature, this model is known as LINTUL3 and simulates crop growth under water-limited and
    nitrogen-limited conditions

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Lintul3.conf"
    __productionlevel__ = "NWLP"
    __cropmodel__ = "LINTUL"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = "CNB"


class FAO_WRSI10_WLP_CWB(Engine):
    """Convenience class for computing actual crop water use using the Water Requirements
    Satisfaction Index with a (modified) FAO WRSI approach.

    see `pcse.engine.Engine` for description of arguments and keywords
    """

    config = "FAO_WRSI.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "FAO_WRSI"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = None


class Lingra10_PP(Engine):
    """Convenience class for running the LINGRA grassland model for potential production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Lingra_PP.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "LINGRA"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = None
    __nitrogenbalance__ = None


class Lingra10_WLP_CWB(Engine):
    """Convenience class for running the LINGRA grassland model for water-limited production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Lingra_WLP_FD.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "LINGRA"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = None


class Lingra10_NWLP_CWB_CNB(Engine):
    """Convenience class for running the LINGRA grassland model for nitrogen and water-limited production.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Lingra_NWLP_FD.conf"
    __productionlevel__ = "NWLP"
    __cropmodel__ = "LINGRA"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = "CNB"


class Wofost81_PP(Engine):
    """Convenience class for running WOFOST8.1 potential production

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_PP.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = None
    __nitrogenbalance__ = None


class Wofost81_WLP_CWB(Engine):
    """Convenience class for running WOFOST8.1  water-limited production using the classic
    waterbalance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_WLP_CWB.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = None


class Wofost81_WLP_MLWB(Engine):
    """Convenience class for running WOFOST8.1  water-limited production using the multi-layer
    waterbalance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_WLP_MLWB.conf"
    __productionlevel__ = "WLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = "MLWB"
    __nitrogenbalance__ = None


class Wofost81_NWLP_CWB_CNB(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the classic waterbalance and classic nitrogen balance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_NWLP_CWB_CNB.conf"
    __productionlevel__ = "NWLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = "CWB"
    __nitrogenbalance__ = "CNB"


class Wofost81_NWLP_MLWB_CNB(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the multi-layer waterbalance and classic nitrogen balance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_NWLP_MLWB_CNB.conf"
    __productionlevel__ = "NWLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = "MLWB"
    __nitrogenbalance__ = "CNB"


class Wofost81_NWLP_MLWB_SNOMIN(Engine):
    """Convenience class for running WOFOST8.1 nutrient and water-limited production
    using the multi-layer waterbalance and the SNOMIN carbon/nitrogen balance.

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Wofost81_NWLP_MLWB_SNOMIN.conf"
    __productionlevel__ = "NWLP"
    __cropmodel__ = "WOFOST"
    __cropmodelversion__ = "8.1"
    __waterbalance__ = "MLWB"
    __nitrogenbalance__ = "SNOMIN"


class Alcepas10_PP(Engine):
    """Convenience class for running the ALCEPAS 1.0 onion model for potential production

    see `pcse.engine.Engine` for description of arguments and keywords
    """
    config = "Alcepas10_PP.conf"
    __productionlevel__ = "PP"
    __cropmodel__ = "ALCEPAS"
    __cropmodelversion__ = "1.0"
    __waterbalance__ = None
    __nitrogenbalance__ = None


# This is to ensure that old code keeps working
Wofost71_PP = Wofost72_PP
Wofost71_WLP_FD = Wofost72_WLP_FD = Wofost72_WLP_CWB
LINTUL3 = Lintul10_NWLP_CWB_CNB
LINGRA_PP = Lingra10_PP
LINGRA_WLP_FD = Lingra10_WLP_CWB
LINGRA_NWLP_FD = Lingra10_NWLP_CWB_CNB
