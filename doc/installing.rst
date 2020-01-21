***************
Installing PCSE
***************

Requirements and dependencies
=============================

PCSE is being developed on Ubuntu Linux 14.04 and Windows 7 using Python 2.7.14 and python 3.6.5.
As Python is a platform independent language, PCSE works equally well on Linux, Windows or Mac OSX.
Before installing PCSE, Python itself must be installed on your system which we will demonstrate
below. PCSE has a number of dependencies on other python packages which are the following:

* Numpy >= 1.6
* SQLalchemy >= 0.8
* PyYAML >= 3.11
* xlrd >= 0.9.3
* xlwt >= 1.0.0
* pandas >= 0.20
* requests >= 2.0.0
* traitlets-pcse == 5.0.0.dev

The last package in the list is a modified version of the `traitlets`_ package which provides some
additional functionality used by PCSE.

.. _Enthought Canopy: https://www.enthought.com/products/canopy/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _PythonXY: https://python-xy.github.io/
.. _HomeBrew: http://brew.sh
.. _traitlets: https://traitlets.readthedocs.io/en/stable/

Setting up your python environment
==================================

A convenient way to set up your python environment for PCSE is through the `Anaconda`_ python distribution.
In the present PCSE Documentation all examples of installing and using PCSE refer to the Windows 7 platform.

First, we suggest you download and install the `MiniConda`_ python distribution which provides a minimum
python environment that we will use to bootstrap a dedicated virtual environment for PCSE. For the rest
of this guide we will assume that you use Windows 7 and install the
32bit miniconda for python 3 (``Miniconda3-latest-Windows-x86.exe``). The virtual environment that
we will create contains not only the dependencies for PCSE, it also includes many other useful packages
such as `IPython`_, `Pandas`_ and the `Jupyter notebook`_. These packages will be used in the Getting Started section
as well.

.. _MiniConda: http://conda.pydata.org/miniconda.html
.. _Pandas: http://pandas.pydata.org/
.. _Jupyter notebook: https://jupyter.org/
.. _IPython: https://ipython.org/

After installing MiniConda you should open a command box and check that conda is installed properly:

.. code-block:: doscon

    C:\>conda info
    Current conda install:

                 platform : win-32
            conda version : 4.0.5
      conda-build version : not installed
           python version : 3.5.1.final.0
         requests version : 2.9.1
         root environment : C:\Miniconda3  (writable)
      default environment : C:\Miniconda3
         envs directories : C:\Miniconda3\envs
            package cache : C:\Miniconda3\pkgs
             channel URLs : https://repo.continuum.io/pkgs/free/win-32/
                            https://repo.continuum.io/pkgs/free/noarch/
                            https://repo.continuum.io/pkgs/pro/win-32/
                            https://repo.continuum.io/pkgs/pro/noarch/
              config file : None
        is foreign system : False

Now we will use a Conda environment file to recreate the python environment that we use to develop and run
PCSE. First you should download the conda environment file which comes in two flavours, an
environment for running PCSE on python 3 (:download:`downloads/py3_pcse.yml`) and one for python 2
(:download:`downloads/py2_pcse.yml`). Both environments include the Jupyter notebook and IPython which are
needed for running the `getting started` section and the example notebooks. Save the environment file
on a temporary location such as ``d:\temp\make_env\``. We will now create a dedicated virtual environment
using the command ``conda env create`` and tell conda to use the environment file for python3 with the
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

Installing and testing PCSE
===========================

The easiest way to install PCSE is through the python package index (`PyPI`_).
Installing from PyPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip` which searches
the python package index for a package, downloads and installs it into your python
environment:

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

If you are wondering what the difference between `pip` and `conda` are than have a look
`here <https://stackoverflow.com/questions/20994716/what-is-the-difference-between-pip-and-conda#20994790>`_

If you want to develop with or contribute to PCSE, than you should fork the `PCSE
repository`_ on GitHub and get a local copy of PCSE using `git clone`. See the help on github_
and for Windows/Mac users the `GitHub Desktop`_ application.

.. _GitHub Desktop: https://desktop.github.com/
.. _GitHub: https://help.github.com/
.. _PCSE repository: https://github.com/ajwdewit/pcse
.. _PyPI: https://pypi.python.org/pypi/PCSE

To guarantee its integrity, the PCSE package includes a number of self
tests that test individual components as well as the entire simulation. These tests
verify that the output produced by the different components matches with the
expected outputs. Test data for the individual components can be found
in the `pcse.tests.test_data` package, while the test data for the entire chain
is stored in an SQLite database (pcse.db). This database can be found under
`.pcse` in your home folder and will be automatically created when importing
PCSE for the first time. When you delete the database file manually it will be
recreated next time you import PCSE.

For testing the PCSE package we need to start python and import pcse:

.. code-block:: doscon

    (py3_pcse) D:\temp\make_env>python
    Python 3.6.5 (default, Aug 14 2018, 19:12:50) [MSC v.1900 32 bit (Intel)] :: Anaconda, Inc. on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db ... OK
    >>>

Next, the tests can be executed by calling the `test()` function at the top of the package::


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
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_lintul3.TestLINTUL3_SpringWheat) ... ok
    runTest (pcse.tests.test_wofost_npk.TestWOFOSTNPK_WinterWheat) ... ok

    ----------------------------------------------------------------------
    Ran 32 tests in 54.306s

    OK
    >>>

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred. Note that the results may deviate from the output above
when tests were added or removed.

Moreover, SQLAlchemy may complain with a warning that can be safely ignored::

    C:\Miniconda3\envs\py3_pcse\lib\site-packages\sqlalchemy\sql\sqltypes.py:603: SAWarning:
    Dialect sqlite+pysqlite does *not* support Decimal objects natively, and SQLAlchemy must
    convert from floating point - rounding errors and other issues may occur. Please consider
    storing Decimal numbers as strings or integers on this platform for lossless storage.

