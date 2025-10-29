Python Crop Simulation Environment - PCSE
=========================================

PCSE is a framework developed for implementing crop simulation models developed in
Wageningen. Many of the Wageningen crop simulation models were originally developed
in FORTRAN77 or using the FORTRAN Simulation Translator (`FST`_). Although this
approach has yielded high quality models with high numerical performance, the
inherent limitations of models written in FORTRAN is also becoming increasingly evident:

* The structure of the models is often rather monolithic and the different parts are
  very tightly coupled. Replacing parts of the model with another simulation approach
  is not easy.
* The models rely on file-based I/O which is difficult to change. For example,
  interfacing with databases is complicated in FORTRAN.
* In general, with low-level languages like FORTRAN, simple things already take many
  lines of code and mistakes are easily made, particularly by agronomists and crop
  scientist that have limited experience in developing or adapting software.

To overcome many of the limitations above, the Python Crop Simulation Environment
(PCSE) was developed which provides an environment for developing simulation models
as well as a number of implementations of crop simulation models. PCSE is written
in pure python code which makes it more flexible, easier to modify and extensible
allowing easy interfacing with databases, graphical user interfaces, visualization
tools and numerical/statistical packages. PCSE has several interesting features:

* Implementation in pure python with dependencies only on popular packages available from
  the Python Package Index (PyPI) (`Pydispather`, `PyYAML`, `pandas`, `Openpyxl` and  `requests`)

* Modular design allowing you to add or change components relatively quickly with
  a simple but powerful approach to communicate variables between modules.

* Similar to `FST`_, it enforces good model design by explicitly separating parameters,
  rate variables and state variables. Moreover PCSE takes care of the module
  initialization, calculation of rates of changes, updating of state variables
  and actions needed to finalize the simulation.

* Input/Output is completely separated from the simulation model itself. Therefore
  PCSE models can easily read from and write to text files, databases and scientific
  formats such as HDF or NetCDF.

* Tools are available for reading parameter and weather files from existing models to
  have as much backward compatibility as possible.

* An `AgroManager` module which allows to define the agromanagement actions that
  happen on a farmers field. Such actions can be specified as events based on
  time or model state.

* Built-in testing of program modules ensuring integrity of the system.

To contribute to PCSE, you can fork your own copy at https://github.com/ajwdewit/pcse

Full documentation is available on https://pcse.readthedocs.io


Testing PCSE
------------

The PCSE package has some built-in tests that can used to test if any PCSE installation is
producing the correct outputs::

    >>> pcse.test()
    runTest (pcse.tests.test_abioticdamage.Test_FROSTOL) ... ok
    runTest (pcse.tests.test_partitioning.Test_DVS_Partitioning) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_PotentialEvapotranspiration) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok

    ...

    runTest (pcse.tests.test_wofost.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost80.TestWOFOST80_Potential_WinterWheat) ... ok
    runTest (pcse.tests.test_wofost80.TestWOFOST80_WaterLimited_WinterWheat) ... ok

    ----------------------------------------------------------------------
    Ran 32 tests in 39.809s

    OK

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred. Note that the results may deviate from the output above
when tests were added or removed.

On top of the built-in tests, a larger suite of tests is available in the
git repository of PCSE. The latter also includes tests of the LINGRA model
which are not included in the internal tests. The tests can be execute through
the `tests` package::

    (py3_pcse) $ python -m tests
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_potentialproduction_wofost72_01.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_potentialproduction_wofost72_11.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_potentialproduction_wofost72_21.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_potentialproduction_wofost72_31.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_potentialproduction_wofost72_41.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_waterlimitedproduction_wofost72_01.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_waterlimitedproduction_wofost72_11.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_waterlimitedproduction_wofost72_21.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_waterlimitedproduction_wofost72_31.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_waterlimitedproduction_wofost72_41.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Belgium-Michamps-1986_PP.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Netherlands-Zegveld-1986_PP.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Belgium-Michamps-1986_WLP.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Netherlands-Zegveld-1986_WLP.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Belgium-Michamps-1986_NWLP.yaml) ... ok
    runTest (tests.run_tests./home/wit015/Sources/python/pcse/tests/test_data/test_LINGRA_Netherlands-Zegveld-1986_NWLP.yaml) ... ok

    ----------------------------------------------------------------------
    Ran 16 tests in 101.956s

    OK

By default this runs a limited selection of tests. The full test suite can be run with::

    (py3_pcse) $ python -m tests --full

But this will take at least 30 minutes to complete.


Comparing PCSE models against experiments
-----------------------------------------

Starting with PCSE 5.5, there is an additional folder `exp` inside the repository which contains experimental
data which can be used to compare the results from a PCSE model against. Experiments are collected in an
'experimental collection' which contains references to experiments that belong together. For example, all
experiments for potato for a given variety. Currently, the available experiments are limited to grassland for the
LINGRA model and consist of two collections. One for grassland under irrigated conditions and one for rain-fed
conditions. Tt is expected that more experimental data will be collected and stored here in order to have a
reference set to compare model results.

Running the experiments is similar to running the unit tests::

    (py3_pcse) $ python -m exp
    Writing expriment results to: /tmp/exp_results
    Processing collection for Rye grass: Potential
      - Processing experiment: LINGRA_FAO/LINGRA_FAO_experiment_000_UK2_1982.yaml
      - Processing experiment: LINGRA_FAO/LINGRA_FAO_experiment_004_SW1_1983.yaml
      - Processing experiment: LINGRA_FAO/LINGRA_FAO_experiment_006_SW1_1984.yaml

    ...

This will generate figures of simulated vs observed data in order to assess how the model
performs against experimental data. In the future, this will be extended to include a
report with error values.


.. _FST: http://models.pps.wur.nl/sites/models.pps.wur.nl/files/FST%203.pdf