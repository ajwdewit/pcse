What's new in PCSE 5.2
======================

PCSE version 5.2 brings the following new features:

- The LINTUL3 model has been implemented in PCSE. LINTUL3 is a simple crop growth model for simulating
  growth conditions under water-limited and nitrogen-limited conditions.
- A new module for N/P/K limitations in WOFOST was implemented allowing to simulate the impact of N/P/K
  limitations on crop growth in WOFOST.
- A new :ref:`AgroManager <Agromanagement>` which greatly enhances the way that AgroManagement can be handled in PCSE.
  The new agromanager
  can elegantly combine cropping calendars, timed events and state events also within rotations over several cropping
  campaigns. The AgroManager uses a new format based on YAML to store agromanagement definitions.

What's new in PCSE 5.1
======================

PCSE version 5.1 brings the following new features:

- Support for reading input data (weather, soil, crop parameters) from a CGMS11 database. CGMS is the acronym for
  Crop Growth Monitoring System and was developed by Alterra in cooperation with the MARS unit of the Joint Research
  Centre for crop monitoring and yield forecasting in Europe. It uses a database structure for storing weather
  data and model simulation results which can be read by PCSE. See the MARSwiki_ for the database definition.
- The ExcelWeatherDataProvider: Before PCSE 5.2 the only file-based format for weather data was the CABO weather format
  read by the :ref:`CABOWeatherDataProvider <CABOWeatherDataProvider>`. Althought the format is well documented,
  creating CABO weather files is a bit cumbersome as for each year a new file has to be created and mistakes are
  easily made. Therefore, the :ref:`ExcelWeatherDataProvider <ExcelWeatherDataProvider>` was created that
  reads its input from a Microsoft Excel file. See here for an example of an Excel weather file: :download:`nl1.xlsx`.


.. _MARSwiki: http://marswiki.jrc.ec.europa.eu/agri4castwiki/index.php/Appendix_5:_CGMS_tables