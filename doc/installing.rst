Installing PCSE
===============

Requirements and dependencies
-----------------------------

PCSE is being developed on Ubuntu Linux 18.04 and Windows 10 using Python 3.7 and Python 3.8
As Python is a platform independent language, PCSE works equally well on Linux, Windows or Mac OSX.
Before installing PCSE, Python itself must be installed on your system which we will demonstrate
below. PCSE has a number of dependencies on other Python packages which are the following::

- SQLAlchemy>=0.8.0
- PyYAML>=3.11
- xlrd>=0.9.3
- openpyxl>=3.0
- requests>=2.0.0
- pandas>=0.20
- traitlets-pcse==5.0.0.dev


The last package in the list is a modified version of the `traitlets`_ package which provides some
additional functionality used by PCSE.

.. _Enthought Canopy: https://www.enthought.com/products/canopy/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _PythonXY: https://python-xy.github.io/
.. _HomeBrew: http://brew.sh
.. _traitlets: https://traitlets.readthedocs.io/en/stable/

Setting up your Python environment
----------------------------------

A convenient way to set up your Python environment for PCSE is through the `Anaconda`_ Python distribution.
In the present PCSE Documentation all examples of installing and using PCSE refer to the Windows 10 platform.

First, we suggest you download and install the `MiniConda`_ Python distribution which provides a minimum
Python environment that we will use to bootstrap a dedicated environment for PCSE. For the rest
of this guide we will assume that you use Windows 10 and install the
64bit miniconda for Python 3 (``Miniconda3-latest-Windows-x86_64.exe``). The environment that
we will create contains not only the dependencies for PCSE, it also includes many other useful packages
such as `IPython`_, `Pandas`_ and the `Jupyter notebook`_. These packages will be used in the Getting Started section
as well.

.. _MiniConda: http://conda.pydata.org/miniconda.html
.. _Pandas: http://pandas.pydata.org/
.. _Jupyter notebook: https://jupyter.org/
.. _IPython: https://ipython.org/

After installing MiniConda you should open a command box and check that conda is installed properly:

.. code-block:: doscon

    (py3_pcse) C:\>conda info

             active environment : py3_pcse
            active env location : C:\data\Miniconda3\envs\py3_pcse
                    shell level : 3
               user config file : C:\Users\wit015\.condarc
         populated config files : C:\Users\wit015\.condarc
                  conda version : 4.9.2
            conda-build version : not installed
                 python version : 3.8.5.final.0
               virtual packages : __win=0=0
                                  __archspec=1=x86_64
               base environment : C:\data\Miniconda3  (writable)
                   channel URLs : https://conda.anaconda.org/conda-forge/win-64
                                  https://conda.anaconda.org/conda-forge/noarch
                                  https://repo.anaconda.com/pkgs/main/win-64
                                  https://repo.anaconda.com/pkgs/main/noarch
                                  https://repo.anaconda.com/pkgs/r/win-64
                                  https://repo.anaconda.com/pkgs/r/noarch
                                  https://repo.anaconda.com/pkgs/msys2/win-64
                                  https://repo.anaconda.com/pkgs/msys2/noarch
                  package cache : C:\data\Miniconda3\pkgs
                                  C:\Users\wit015\.conda\pkgs
                                  C:\Users\wit015\AppData\Local\conda\conda\pkgs
               envs directories : C:\data\Miniconda3\envs
                                  C:\Users\wit015\.conda\envs
                                  C:\Users\wit015\AppData\Local\conda\conda\envs
                       platform : win-64
                     user-agent : conda/4.9.2 requests/2.24.0 CPython/3.8.5 Windows/10 Windows/10.0.18362
                  administrator : False
                     netrc file : None
                   offline mode : False

Now we will use a Conda environment file to recreate the Python environment that we use to develop and run
PCSE. First you should download the conda environment file which comes in two flavours, an
environment for running PCSE on Python 3 (:download:`downloads/py3_pcse.yml`) and one for Python 2
(:download:`downloads/py2_pcse.yml`). It is strongly recommended to use the Python 3 version as Python 2
is not maintained anymore. Both environments include the Jupyter notebook and IPython which are
needed for running the `getting started` section and the example notebooks. Save the environment file
on a temporary location such as ``d:\temp\make_env\``. We will now create a dedicated virtual environment
using the command ``conda env create`` and tell conda to use the environment file for Python 3 with the
option ``-f p3_pcse.yml`` as show below:

.. code-block:: doscon

    (C:\Miniconda3) D:\temp\make_env>conda env create -f py3_pcse.yml
    Fetching package metadata .............
    Solving package specifications: .
    intel-openmp-2 100% |###############################| Time: 0:00:00   6.39 MB/s

    ... Lots of output here

    Installing collected packages: traitlets-pcse
    Successfully installed traitlets-pcse-5.0.0.dev0
    #
    # To activate this environment, use:
    # > activate py3_pcse
    #
    # To deactivate an active environment, use:
    # > deactivate
    #
    # * for power-users using bash, you must source
    #

You can then activate your environment (note the addition of ``(py3_pcse)`` on your command prompt):

.. code-block:: doscon

    D:\temp\make_env>activate py3_pcse
    Deactivating environment "C:\Miniconda3"...
    Activating environment "C:\Miniconda3\envs\py3_pcse"...

    (py3_pcse) D:\temp\make_env>

Installing PCSE
---------------

The easiest way to install PCSE is through the Python package index (`PyPI`_).
Installing from PyPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip` which searches
the Python package index for a package, downloads and installs it into your Python
environment (example below for PCSE 5.4):

.. code-block:: doscon

    (py3_pcse) D:\temp\make_env>pip install pcse

    Collecting pcse
      Downloading https://files.pythonhosted.org/packages/8c/92/d4444cce1c58e5a96f4d6dc9c0e042722f2136df24a2750352e7eb4ab053/PCSE-5.4.0.tar.gz (791kB)
        100% |¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦| 798kB 1.6MB/s
    Requirement already satisfied: numpy>=1.6.0 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (1.15.1)
    Requirement already satisfied: SQLAlchemy>=0.8.0 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (1.2.11)
    Requirement already satisfied: PyYAML>=3.11 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (3.13)
    Requirement already satisfied: xlrd>=0.9.3 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (1.1.0)
    Requirement already satisfied: xlwt>=1.0.0 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (1.3.0)
    Requirement already satisfied: requests>=2.0.0 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (2.19.1)
    Requirement already satisfied: pandas>=0.20 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (0.23.4)
    Requirement already satisfied: traitlets-pcse==5.0.0.dev in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pcse) (5.0.0.dev0)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from requests>=2.0.0->pcse) (3.0.4)
    Requirement already satisfied: idna<2.8,>=2.5 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from requests>=2.0.0->pcse) (2.7)
    Requirement already satisfied: certifi>=2017.4.17 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from requests>=2.0.0->pcse) (2018.8.24)
    Requirement already satisfied: urllib3<1.24,>=1.21.1 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from requests>=2.0.0->pcse) (1.23)
    Requirement already satisfied: python-dateutil>=2.5.0 in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pandas>=0.20->pcse) (2.7.3)
    Requirement already satisfied: pytz>=2011k in c:\miniconda3\envs\py3_pcse\lib\site-packages (from pandas>=0.20->pcse) (2018.5)
    Requirement already satisfied: six in c:\miniconda3\envs\py3_pcse\lib\site-packages (from traitlets-pcse==5.0.0.dev->pcse) (1.11.0)
    Requirement already satisfied: decorator in c:\miniconda3\envs\py3_pcse\lib\site-packages (from traitlets-pcse==5.0.0.dev->pcse) (4.3.0)
    Requirement already satisfied: ipython-genutils in c:\miniconda3\envs\py3_pcse\lib\site-packages (from traitlets-pcse==5.0.0.dev->pcse) (0.2.0)
    Building wheels for collected packages: pcse
      Running setup.py bdist_wheel for pcse ... done
      Stored in directory: C:\Users\wit015\AppData\Local\pip\Cache\wheels\2f\e6\2c\3952ff951dffea5ab2483892edcb7f9310faa319d050d3be6c
    Successfully built pcse
    twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.
    mkl-random 1.0.1 requires cython, which is not installed.
    mkl-fft 1.0.4 requires cython, which is not installed.
    Installing collected packages: pcse
    Successfully installed pcse-5.4.0

If you are wondering what the differences between `pip` and `conda` are, then have a look
`here <https://stackoverflow.com/questions/20994716/what-is-the-difference-between-pip-and-conda#20994790>`_

If you want to develop with or contribute to PCSE, then you should fork the `PCSE
repository`_ on GitHub and get a local copy of PCSE using `git clone`. See the help on GitHub_
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

For running the internal tests of the PCSE package we need to start Python and import pcse:

.. code-block:: doscon

    (py3_pcse) D:\temp\make_env>python
    Python 3.6.5 (default, Aug 14 2018, 19:12:50) [MSC v.1900 32 bit (Intel)] :: Anaconda, Inc. on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db ... OK
    >>>

Next, the tests can be executed by calling the `test()` function at the top of the package::

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
    runTest (pcse.tests.test_wofost.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok
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

Moreover, SQLAlchemy may complain with a warning that can be safely ignored::

    C:\Miniconda3\envs\py3_pcse\lib\site-packages\sqlalchemy\sql\sqltypes.py:603: SAWarning:
    Dialect sqlite+pysqlite does *not* support Decimal objects natively, and SQLAlchemy must
    convert from floating point - rounding errors and other issues may occur. Please consider
    storing Decimal numbers as strings or integers on this platform for lossless storage.
