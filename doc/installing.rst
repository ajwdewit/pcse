***************
Installing PCSE
***************

Requirements and dependencies
=============================

PCSE is being developed on Ubuntu Linux 14.04 and Windows 7 using Python 2.7.11 and is known to work with
the 3.x series (using the 2to3 tool). As Python is a platform independent language, PCSE
works equally well on Linux, Windows or Mac OSX.
Before installing PCSE, Python itself must be installed on your system. Most Linux systems provide
Python through the native package manager. For Windows users the most straightforward approach for installing
Python is through one of the prepackaged Python distributions such as `Enthought Canopy`_,
`Anaconda`_ or `PythonXY`_. The advantage of the prepackaged distributions is that they provide a working
version of Numpy out-of-the-box which can be difficult to install on Windows. Mac OSX users can most easily
install Python and Numpy using `HomeBrew`_, i.e brew install Python Numpy

The dependencies of PCSE are the following:

* Numpy >= 1.6
* SQLalchemy >= 0.8
* PyYAML >= 3.11
* tabulate >= 0.7.5
* xlrd >= 0.9.3
* xlwt >= 1.0.0

.. _Enthought Canopy: https://www.enthought.com/products/canopy/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _PythonXY: https://python-xy.github.io/
.. _HomeBrew: http://brew.sh

Setting up your python environment
==================================

A convenient way to set up your python environment for PCSE is through the `Anaconda`_ python distribution.
In the present PCSE Documentation all examples of installing and using PCSE refer to Windows 7 platform.

First, we suggest you download and install the `MiniConda`_ python distribution which provides a minimum
python environment that we will use to bootstrap a dedicated virtual environment for PCSE. For the rest
of this guide we will assume that you use Windows 7 and install the
32bit miniconda for python 3 (``Miniconda3-latest-Windows-x86.exe``). The virtual enviroment that
we will create contains not only the dependencies for PCSE, it also includes many other useful packages
such as `IPython`_ `Pandas`_ and the `Jupyter notebook`_. These packages will be used in the Getting Started section
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
PCSE. First you should download the conda environment file :download:`downloads/py2_pcse.yml` and save
it on a temporary location such as ``d:\temp\make_env\``. We will now create a dedicated virtual environment
using the command ``conda env create`` and tell conda to use the environment file with the option ``-f p2_pcse.yml``
as show below:

.. code-block:: doscon

    D:\temp\make_env>conda env create -f py2_pcse.yml
    Fetching package metadata: ....
    .Solving package specifications: .........
    Linking packages ...
    [      COMPLETE      ]|##################################################| 100%
    Requirement already satisfied (use --upgrade to upgrade): backports.ssl-match-hostname==3.4.0.2 in c:\miniconda3\envs\py2_pcse_b\lib\site-packages
    Requirement already satisfied (use --upgrade to upgrade): line-profiler==1.0 in c:\miniconda3\envs\py2_pcse_b\lib\site-packages
    Requirement already satisfied (use --upgrade to upgrade): sphinx-rtd-theme==0.1.7 in c:\miniconda3\envs\py2_pcse_b\lib\site-packages
    Collecting tabulate==0.7.5
    Installing collected packages: tabulate
    Successfully installed tabulate-0.7.5
    #
    # To activate this environment, use:
    # > activate py2_pcse
    #

You can then activate your environment (note the addition of ``[py2_pcse]`` on your command prompt) and
run an upgrade to upgrade all packages to the latest versions:

.. code-block:: doscon

    D:\temp\make_env>activate py2_pcse
    Deactivating environment "C:\Miniconda3"...
    Activating environment "C:\Miniconda3\envs\py2_pcse"...

    [py2_pcse] D:\temp\make_env>conda upgrade --all
    Fetching package metadata: ....
    Solving package specifications: .........

    Package plan for installation in environment C:\Miniconda3\envs\py2_pcse:

    The following packages will be downloaded:

        package                    |            build
        ---------------------------|-----------------
        vs2008_runtime-9.00.30729.1|                1         1.1 MB
        zlib-1.2.8                 |            vc9_3          95 KB

    ... Lots of output here

    Proceed ([y]/n)? y

    Fetching packages ...
    vs2008_runtime 100% |###############################| Time: 0:00:02 518.43 kB/s

    ... Lots of output here

    pyqt-4.11.4-py 100% |###############################| Time: 0:00:04 571.51 kB/s
    tornado-4.3-py 100% |###############################| Time: 0:00:01 387.81 kB/s
    Extracting packages ...
    [      COMPLETE      ]|##################################################| 100%
    Unlinking packages ...
    [      COMPLETE      ]|##################################################| 100%
    Linking packages ...
            1 file(s) copied.######################################          |  81%
    [      COMPLETE      ]|##################################################| 100%

    [py2_pcse] D:\temp\make_env>


Installing and testing PCSE
===========================

The easiest way to install PCSE is through the python package index (`PyPI`_).
Installing from PyPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip` which searches
the python package index for a package, downloads and installs it into your python
environment:

.. code-block:: doscon

    [py2_pcse] D:\temp\make_env>pip install PCSE
    Collecting PCSE
    Requirement already satisfied (use --upgrade to upgrade): numpy>=1.6.0 in c:\miniconda3\envs\py2_pcse\lib\site-packages (from PCSE)
    Requirement already satisfied (use --upgrade to upgrade): xlrd>0.9.0 in c:\miniconda3\envs\py2_pcse\lib\site-packages (from PCSE)
    Requirement already satisfied (use --upgrade to upgrade): tabulate>=0.7.0 in c:\miniconda3\envs\py2_pcse\lib\site-packages (from PCSE)
    Requirement already satisfied (use --upgrade to upgrade): SQLAlchemy>=0.8.0 in c:\miniconda3\envs\py2_pcse\lib\site-packages (from PCSE)
    Installing collected packages: PCSE
    Successfully installed PCSE-5.2

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

    [py2_pcse_b] D:\temp\make_env>python
    Python 2.7.11 |Continuum Analytics, Inc.| (default, Mar  4 2016, 15:18:41) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://anaconda.org
    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db
    >>>


Next, the tests can be executed by calling the `test()` function at the top of the package::

    >>> pcse.test()
    runTest (pcse.tests.test_abioticdamage.Test_FROSTOL) ... ok
    runTest (pcse.tests.test_assimilation.Test_WOFOST_Assimilation) ... ok
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
    runTest (pcse.tests.test_wofost.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_lintul3.TestLINTUL3_SpringWheat) ... ok
    runTest (pcse.tests.test_wofost_npk.TestWOFOSTNPK_WinterWheat) ... ok

    ----------------------------------------------------------------------
    Ran 33 tests in 57.472s

    OK

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred. Note that the results may deviate from the output above
because one or more tests may have been temporarily disabled (skipped) often
due to problems with the test.

Moreover, SQLAlchemy may complain with a warning that can be safely ignored::

     /usr/lib/python2.7/dist-packages/sqlalchemy/types.py:307: SAWarning:
     Dialect sqlite+pysqlite does *not* support Decimal objects natively, and
     SQLAlchemy must convert from floating point - rounding errors and other
     issues may occur. Please consider storing Decimal numbers as strings or
     integers on this platform for lossless storage.
         d[coltype] = rp = d['impl'].result_processor(dialect, coltype)

