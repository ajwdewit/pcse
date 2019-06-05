What's new in PCSE 5.4
======================

PCSE 5.4 has the following new features:

- PCSE is now fully compatible with python3 (>3.4) while still remaining compatibility with python 2.7.14
- The NASAPOWERWeatherDataProvider has been upgraded to take the new API into account



What's new in PCSE 5.3
======================

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

What's new in PCSE 5.2
======================

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

What's new in PCSE 5.1
======================

PCSE version 5.1 brings the following new features:

- Support for reading input data (weather, soil, crop parameters) from a CGMS12 database. CGMS is the acronym for
  Crop Growth Monitoring System and was developed by Alterra in cooperation with the MARS unit of the Joint Research
  Centre for crop monitoring and yield forecasting in Europe. It uses a database structure for storing weather
  data and model simulation results which can be read by PCSE. See the MARSwiki_ for the database definition.
- The ExcelWeatherDataProvider: Before PCSE 5.2 the only file-based format for weather data was the CABO weather format
  read by the :ref:`CABOWeatherDataProvider <CABOWeatherDataProvider>`. Althought the format is well documented,
  creating CABO weather files is a bit cumbersome as for each year a new file has to be created and mistakes are
  easily made. Therefore, the :ref:`ExcelWeatherDataProvider <ExcelWeatherDataProvider>` was created that
  reads its input from a Microsoft Excel file. See here for an example of an Excel weather file: :download:`downloads/nl1.xlsx`.


.. _MARSwiki: http://marswiki.jrc.ec.europa.eu/agri4castwiki/index.php/Appendix_5:_CGMS_tables