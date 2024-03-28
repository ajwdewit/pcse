Models available in PCSE
========================

The following table lists the models that are available in PCSE and can be imported from `pcse.models` package.
For each model it lists the production level and some of the features included in the models. The model name for
most models is made up from a set of codes which follow: <modelname><version>_<productionlevel>_<waterbalance>_<nitrogenbalance>

.. Table generated using: https://tableconvert.com/restructuredtext-generator

=========================== ============================ ============ ====================== ============ =============== ===========
 Model                       Production level             CO2 impact   Biomass reallocation   N dynamics   Water balance   N balance
=========================== ============================ ============ ====================== ============ =============== ===========
 Wofost72_Pheno              Phenology only                                                                N/A             N/A
 Wofost72_PP                 Potential                                                                     N/A             N/A
 Wofost72_WLP_CWB            Water-limited                                                                 Classic         N/A
 Wofost73_PP                 Potential                    X            X                                   N/A             N/A
 Wofost73_WLP_CWB            Water-limited                X            X                                   Classic         N/A
 Wofost73_WLP_MLWB           Water-limited                X            X                                   Multi-layer     N/A
 Wofost81_PP                 Potential                    X            X                      X            N/A             N/A
 Wofost81_WLP_CWB            Water-limited                X            X                      X            Classic         N/A
 Wofost81_WLP_MLWB           Water-limited                X            X                      X            Multi-layer     N/A
 Wofost81_NWLP_CWB_CNB       Water and Nitrogen limited   X            X                      X            Classic         Classic
 Wofost81_NWLP_MLWB_CNB      Water and Nitrogen limited   X            X                      X            Multi-layer     Classic
 Wofost81_NWLP_MLWB_SNOMIN   Water and Nitrogen limited   X            X                      X            Multi-layer     SNOMIN
 LINGRA_PP                   Potential                    X                                                N/A             N/A
 LINGRA_WLP_CWB              Water-limited                X                                                Classic         N/A
 LINGRA_NWLP_CWB_CNB         Water and Nitrogen limited   X                                   X            Classic         Classic
 LINTUL3_NWLP_CWB_CNB        Water and Nitrogen limited                                       X            Classic         Classic
 ALCEPAS_PP                  Potential                                                                     N/A             N/A
 FAO_WRSI                    Water-limited                                                                 Classic         N/A
=========================== ============================ ============ ====================== ============ =============== ===========


