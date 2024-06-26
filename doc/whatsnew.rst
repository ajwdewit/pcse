.. include:: abbreviations.txt

#####################################
An overview of new features and fixes
#####################################

**********************
What's new in PCSE 6.0
**********************

PCSE 6.0 has a number of breaking changes compared to the 5.5 release.

First of all, WOFOST 7.3 and 8.1 are released with PCSE 6.0. WOFOST 7.3 includes the atmospheric |CO2| response
and biomass reallocation. WOFOST 8.1 provides the full crop N dynamics. Moreover, this release also includes a
new multi-layer waterbalance and a new Carbon/Nitrogen balance called `SNOMIN` (Soil Nitrogen module for Mineral and
Inorganic Nitrogen). Furthermore, WOFOST 8.0-beta has been deprecated in the PCSE 6.0 release. This also implies dropping
support for the simulation of P/K limitations on crop growth.

WOFOST 7.3 and WOFOST 8.1 are coupled to the new SNOMIN carbon/nitrogen balance and the new multi-layer
waterbalance. Nevertheless the existing (simpler) water and nitrogen balances are also still available and in many
cases highly relevant. For example, for educational purposes or in areas with little data available, a simple approach
may be preferred. Therefore, with WOFOST 7.3 and 8.1 many combinations of crop/soil models are
possible. In order to accommodate these combinations, the models in PCSE have been coded according
to a certain system which follows the following rule::

    <modelname><version>_<productionlevel>_<waterbalance>_<nutrientbalance>

The placeholders can be one of the following:
 - <modelname>: Wofost, Lintul, Lingra or any other model included
 - <version>: the respective model version consisting of a major and minor number, e.g. "72" for Wofost 7.2
 - <productionlevel>: an indication of the agro-ecological production level which can be one of
   "PP": potential production, "WLP": water-limited production or "NWLP": nitrogen and water-limited production
   level.
 - <waterbalance>: the type of waterbalance used which can be either one of "CWB": for classic (simpel) water balance
   or "MLWB" for the more complex multi-layer water balance.
 - <nutrientbalance>: the type of nutrient balance which can be either on of "CNB" for the classic (simpel) nitrogen
   balance or "SNOMIN" for the complex carbon/nitrogen balance.

Depending on the production level, the codes for the water balance and nitrogen balance can be omitted. For example,
WOFOST 7.3 potential production will be coded as "Wofost73_PP", while for the most complex WOFOST 8.1 variant the
model coding will be: "Wofost81_NWLP_MLWB_SNOMIN". The latter combines the WOFOST 8.1 crop dynamics with the
multi-layer water balance and the SNOMIN Carbon/Nitrogen balance.

Further changes to the system consist of:

 - all data providers have been moved to `pcse.input`, except for the CGMS dataproviders residing in `pcse.db`.
 - All WOFOST 7.2 models have been renamed according to the coding schema above. So "pcse.models.Wofost72_WLP_FD" is
   now "pcse.models.Wofost72_WLP_CWB" although the old name will keep working in order not to break old code.



**********************
What's new in PCSE 5.5
**********************

PCSE 5.5 has the following new features:

- WOFOST version 8.0 (beta) has been included which has variants for potential (PP), water-limited (WLP) and
  nutrient + water-limited (NWLP) production. Note that dynamics for N/P/K are included in all model variants
  but for the PP and WLP variants the supply of N/P/K is assumed to be unlimited. Note that this a beta version
  because testing of the N/P/K limited growth against experimental data has so far been limited. Nevertheless,
  the dynamics for N/P/K are based on well known principles from other models and rely on the concept of dilution
  curves that define the maximum, critical and residual N/P/K concentration in the crop.
- A full implementation of the LINGRA and LINGRA-N grassland simulation models are now included. This model allows
  to make estimates of productivity of rye grass.
- WOFOST 7.1 has been upgraded to 7.2, this is mainly to be consistent with the updated system description
  for WOFOST at https://wofost.readthedocs.io. Old code that relies on importing WOFOST 7.1 will keep working
  though.
- The WOFOST 7.2 phenology module can now be imported as a standalone model. This is useful when calibration is
  limited to phenology as it greatly increases the model performance.
- The FAO Water Requirement Satisfaction Index is included as a model.


**********************
What's new in PCSE 5.4
**********************

PCSE 5.4 has the following new features:

- PCSE is now fully compatible with python3 (>3.4) while still remaining compatibility with python 2.7.14
- The NASAPOWERWeatherDataProvider has been upgraded to take the new API into account


**********************
What's new in PCSE 5.3
**********************

PCSE 5.3 has the following new features:

- The WOFOST crop parameters have been reorganized into a new data structure and file format (e.g. YAML)
  and are available from github_. PCSE 5.3 provides the :ref:`YAMLCropDataProvider <YAMLCropDataProvider>`
  to read the new parameters files. The YAMLCropDataProvider works together with the AgroManager for
  specifying parameter sets for crop rotations.
- A new :ref:`CGMSEngine <Engine and models>` that mimics the behaviour of the classic CGMS. This means
  the engine can be run up till a specified date. When maturity or harvest is reached, the value of  all
  state variables will be retained and kept constant until the specified date is reached.
- Caching was added to the CGMS weather data providers, this is particularly useful for repeated
  runs as the weather data only have to be retrieved once from the CGMS database.

Some bugs have been fixed:

- The NASA POWER database moved from http:// to https:// so an update of the NASAPowerWeatherDataProvider
  was needed.
- When running crop rotations it was found that python did not garbage collect the crop simulation objects
  quick enough. This is now fixed with an explicit call to the garbage collector.

.. _github: https://github.com/ajwdewit/WOFOST_crop_parameters

**********************
What's new in PCSE 5.2
**********************

PCSE version 5.2 brings the following new features:

- The LINTUL3 model has been implemented in PCSE. LINTUL3 is a simple crop growth model for simulating
  growth conditions under water-limited and nitrogen-limited conditions.
- A new module for N/P/K limitations in WOFOST was implemented allowing to simulate the impact of N/P/K
  limitations on crop growth in WOFOST.
- A new :ref:`AgroManager <refguide_agromanagement>` which greatly enhances the way that AgroManagement can be handled in PCSE.
  The new agromanager
  can elegantly combine cropping calendars, timed events and state events also within rotations over several cropping
  campaigns. The AgroManager uses a new format based on YAML to store agromanagement definitions.
- The water-limited production simulation with WOFOST now supports irrigation using the new AgroManager.
  An example notebook has been added to explain the different irrigation options.
- Support for reading input data from a CGMS8 and CGMS14 database

Changes in 5.2.5:

- Bug fixes in agromanager causing problems with crop_end_type="earliest" or "harvest"
- Caching was added to the CGMS weather data providers
- Added CGMSEngine that mimics behaviour of the classic CGMS: after the cropping season is over, a call
  to _run() will increase the DAY, but the internal state variables do not change anymore, although they
  are kept available and can be queried and stored in OUTPUT.

**********************
What's new in PCSE 5.1
**********************

PCSE version 5.1 brings the following new features:

- Support for reading input data (weather, soil, crop parameters) from a CGMS12 database. CGMS is the acronym for
  Crop Growth Monitoring System and was developed by WEnR in cooperation with the MARS unit of the Joint Research
  Centre for crop monitoring and yield forecasting in Europe. It uses a database structure for storing weather
  data and model simulation results which can be read by PCSE. See the MARSwiki_ for the database definition.
- The ExcelWeatherDataProvider: Before PCSE 5.2 the only file-based format for weather data was the CABO weather format
  read by the :ref:`CABOWeatherDataProvider <CABOWeatherDataProvider>`. Althought the format is well documented,
  creating CABO weather files is a bit cumbersome as for each year a new file has to be created and mistakes are
  easily made. Therefore, the :ref:`ExcelWeatherDataProvider <ExcelWeatherDataProvider>` was created that
  reads its input from a Microsoft Excel file. See here for an example of an Excel weather file: :download:`downloads/nl1.xlsx`.


.. _MARSwiki: http://marswiki.jrc.ec.europa.eu/agri4castwiki/index.php/Appendix_5:_CGMS_tables