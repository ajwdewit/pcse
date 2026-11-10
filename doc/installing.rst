Installing PCSE
===============

Requirements and dependencies
-----------------------------

PCSE is being developed on Ubuntu Linux 26.04 and Windows 11 using python 3.12 and python 3.14
As Python is a platform independent language, PCSE works equally well on Linux, Windows or Mac OSX.
Before installing PCSE, Python itself must be installed on your system which we will demonstrate
below. PCSE has a number of dependencies on other python packages which are the following::

- PyYAML>=5.1
- openpyxl>=3.0
- requests>=2.0
- pandas>=2.0
- dotmap>=1.3
- pydispatcher>=2.0.7
- traitlets-pcse==5.0.0.dev


The last package in the list is a modified version of the `traitlets`_ package which provides some
additional functionality used by PCSE.

.. _Enthought Canopy: https://www.enthought.com/products/canopy/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _PythonXY: https://python-xy.github.io/
.. _HomeBrew: http://brew.sh
.. _traitlets: https://traitlets.readthedocs.io/en/stable/

Setting up your python environment
----------------------------------

A convenient way to set up your python environment for PCSE is through the `Anaconda`_ python distribution.
In the present PCSE Documentation all examples of installing and using PCSE refer to the Windows 11 platform.

First, we suggest you download and install the `MiniConda`_ python distribution which provides a minimum
python environment that we will use to bootstrap a dedicated environment for PCSE. For the rest
of this guide we will assume that you use Windows 11 and install the
64bit miniconda for python 3 (``Miniconda3-latest-Windows-x86_64.exe``). The environment that
we will create contains not only the dependencies for PCSE, it also includes many other useful packages
such as `IPython`_, `Pandas`_ and the `Jupyter notebook`_. These packages will be used in the Getting Started section
as well.

.. _MiniConda: http://conda.pydata.org/miniconda.html
.. _Pandas: http://pandas.pydata.org/
.. _Jupyter notebook: https://jupyter.org/
.. _IPython: https://ipython.org/

After installing MiniConda you should open a command box and check that conda is installed properly:

.. code-block:: doscon

        (base) C:\>conda info

        active environment : base
        active env location : C:\Data\miniconda3
                shell level : 1
           user config file : C:\Users\wit015\.condarc
     populated config files : C:\data\miniconda3\.condarc
                              C:\data\miniconda3\condarc.d\anaconda-auth.yml
                              C:\Users\wit015\.condarc
              conda version : 26.7.1
        conda-build version : not installed
             python version : 3.14.7.final.0
                     solver : libmamba (default)
           virtual packages : __archspec=1=x86_64_v3
                              __conda=26.7.1=0
                              __win=10.0.26200=0
           base environment : C:\data\miniconda3  (writable)
          conda av data dir : C:\data\miniconda3\etc\conda
      conda av metadata url : None
               channel URLs : https://repo.anaconda.com/pkgs/main/win-64
                              https://repo.anaconda.com/pkgs/main/noarch
                              https://repo.anaconda.com/pkgs/r/win-64
                              https://repo.anaconda.com/pkgs/r/noarch
                              https://repo.anaconda.com/pkgs/msys2/win-64
                              https://repo.anaconda.com/pkgs/msys2/noarch
              package cache : C:\data\nobackup\conda_pkgs
              notices cache : C:\Users\wit015\AppData\Local\conda\conda\Cache\notices
           envs directories : C:\data\miniconda3\envs
                              C:\Users\wit015\.conda\envs
                              C:\Users\wit015\AppData\Local\conda\conda\envs
        temporary directory : C:\Users\wit015\AppData\Local\Temp
                   platform : win-64
                 user-agent : conda/26.7.1 requests/2.34.2 CPython/3.14.7 Windows/11 Windows/10.0.26200 solver/libmamba conda-libmamba-solver/26.7.0 libmambapy/2.3.2 aau/0.8.1 c/. s/. e/. i/.
              administrator : False
                 netrc file : None
               offline mode : False


Now we will use create a separate Conda environment that we use to develop and run PCSE as show below:

.. code-block:: doscon

    (base) C:\>conda create -n py314_pcse python=3.14

        3 channel Terms of Service accepted
        Retrieving notices: done
        Channels:
         - defaults
        Platform: win-64
        Collecting package metadata (repodata.json): done
        Solving environment: done

        ... lots of output here

        Downloading and Extracting Packages:

        Preparing transaction: done
        Verifying transaction: done
        Executing transaction: done
        #
        # To activate this environment, use
        #
        #     $ conda activate py314_pcse
        #
        # To deactivate an active environment, use
        #
        #     $ conda deactivate


You can then activate your environment (note the addition of ``(py314_pcse)`` on your command prompt):

.. code-block:: doscon

    C:>conda activate py314_pcse
    Deactivating environment "C:\Miniconda3"...
    Activating environment "C:\Miniconda3\envs\py314_pcse"...

    (py314_pcse) c:\>

Installing PCSE
---------------

The easiest way to install PCSE is through the python package index (`PyPI`_).
Installing from PyPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip` which searches
the python package index for a package, downloads and installs it into your python
environment (example below for PCSE 6.0.13):

.. code-block:: doscon

    (py314_pcse) C:\>pip install pcse

        Collecting pcse
          Using cached pcse-6.0.13-py3-none-any.whl.metadata (9.6 kB)
        Collecting pydispatcher>=2.0.7 (from pcse)
          Using cached PyDispatcher-2.0.7-py3-none-any.whl.metadata (2.4 kB)
        Collecting PyYAML>=5.1 (from pcse)
          Using cached pyyaml-6.0.3-cp314-cp314-win_amd64.whl.metadata (2.4 kB)

        ... lots of output here

        Installing collected packages: pydispatcher, ipython_genutils, dotmap, urllib3, tzdata, six, PyYAML, numpy, idna, et-xmlfile, decorator, charset_normalizer, certifi, traitlets-pcse, requests, python-dateutil, openpyxl, pandas, pcse
        Successfully installed PyYAML-6.0.3 certifi-2026.10.23 charset_normalizer-3.5.2 decorator-5.3.1 dotmap-1.3.30 et-xmlfile-2.0.0 idna-3.20 ipython_genutils-0.2.0 numpy-2.5.4 openpyxl-3.1.5 pandas-3.1.1 pcse-6.0.13 pydispatcher-2.0.7 python-dateutil-4.0.0 requests-2.34.2 six-1.17.0 traitlets-pcse-5.0.0.dev0 tzdata-2026.6 urllib3-2.8.0

        (py314_pcse) C:\>


If you want to develop with or contribute to PCSE, than you should fork the `PCSE
repository`_ on GitHub and get a local copy of PCSE using `git clone`. See the help on github_
and for Windows/Mac users the `GitHub Desktop`_ application.

.. _GitHub Desktop: https://desktop.github.com/
.. _GitHub: https://help.github.com/
.. _PCSE repository: https://github.com/ajwdewit/pcse
.. _PyPI: https://pypi.python.org/pypi/PCSE

Testing PCSE
------------

To guarantee its integrity, the PCSE package includes a limited number of internal
tests that are installed automatically with PCSE. In addition, the PCSE
git repository has a large number of the tests in the `test` folder which do a more
thorough job in testing but will take a long time to complete (e.g. an hour or more).
The internal tests present users with a quick way to ensure that the output produced
by the different components matches with the expected outputs. While the full test
suite is useful for developers only.

Test data for the internal tests can be found in the `pcse.tests.test_data` package as
well as in an SQLite database (pcse.db). This database can be found under
`.pcse` in your home folder and will be automatically created when importing
PCSE for the first time. When you delete the database file manually it will be
recreated next time you import PCSE.

For running the internal tests of the PCSE package we need to start python and import pcse:

.. code-block:: doscon

    (py314_pcse) C:\>python
    Python 3.14.10 | packaged by Anaconda, Inc. | (main, Oct 27 2026, 01:09:10) [MSC v.1942 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db ... OK
    >>>

Next, the tests can be executed by calling the `test()` function at the top of the package:

.. code-block:: doscon

    >>> pcse.test()
    runTest (pcse.tests.test_abioticdamage.Test_FROSTOL) ... ok
    runTest (pcse.tests.test_partitioning.Test_DVS_Partitioning) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_PotentialEvapotranspiration) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_WaterLimitedEvapotranspiration1) ... ok
    runTest (pcse.tests.test_evapotranspiration.Test_WaterLimitedEvapotranspiration2) ... ok
    runTest (pcse.tests.test_respiration.Test_WOFOSTMaintenanceRespiration) ... ok
    runTest (pcse.tests.test_penmanmonteith.Test_PenmanMonteith1) ... ok
    runTest (pcse.tests.test_penmanmonteith.Test_PenmanMonteith2) ... ok
    runTest (pcse.tests.test_penmanmonteith.Test_PenmanMonteith3) ... ok
    runTest (pcse.tests.test_penmanmonteith.Test_PenmanMonteith4) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager1) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager2) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager3) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager4) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager5) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager6) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager7) ... ok
    runTest (pcse.tests.test_agromanager.TestAgroManager8) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedWinterWheat) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialWinterWheat) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost72.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost72.TestWaterlimitedSpringBarley) ... ok

    ----------------------------------------------------------------------
    Ran 30 tests in 22.482s

    OK

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred. Note that the results may deviate from the output above
when tests were added or removed.
