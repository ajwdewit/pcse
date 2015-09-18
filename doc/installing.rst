***************
Installing PCSE
***************

Requirements and dependencies
=============================

PCSE is being developed on Ubuntu Linux 14.04 and Windows 7 using Python 2.7.9 and is known to work with
the 3.x series (using the 2to3 tool). As Python is a platform independent language, PCSE
works equally well on Linux, Windows or Mac OSX.
Before installing PCSE, Python itself must be installed on your system. Most Linux systems provide
Python through the native package manager. For Windows users the most straightforward approach for installing
Python is through one of the prepackaged Python distributions such as `Enthought Canopy`_,
`Anaconda`_ or `PythonXY`_. The advantage of the prepackaged distributions is that they provide a working
version of Numpy out-of-the-box which can be difficult to install on Windows. Mac OSX users can most easily
install Python and numpy using `HomeBrew`_, i.e brew install python numpy

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

How to install
==============

PCSE can be installed in many different ways and the best depends on your
requirements. The most convenient option to install PCSE is through the Python Package
Index (PyPI). Installing from PyPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip`, but this
will only work when you have write access on the Python site-packages
folder. If you do not have this permission you should install PCSE into a
`virtual environment`_.::

    $ pip install PCSE
    Downloading/unpacking PCSE
      Running setup.py install for PCSE
      ...
      ...
    Successfully installed PCSE
    Cleaning up...
    $

.. _virtual environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/

The second option is useful if you want to develop or contribute to PCSE.
In that case you should fork my `PCSE
repository`_ on GitHub and get a local copy of PCSE using `git clone`. See the help on github_
and for Windows/Mac users the `GitHub Desktop`_ application.

.. _GitHub Desktop: https://desktop.github.com/
.. _GitHub: https://help.github.com/
.. _PCSE repository: https://github.com/ajwdewit/pcse

Finally, the last option is to download the PCSE package as a zip file from GitHub
using the link `here`_. Just unzip the package at a suitable location.
Note that the top directory in the zip file is `pcse-<branchname>`.
The actual PCSE is package is inside this folder and needs to be put on your file system.
For running the test included with PCSE, we assume that PCSE was downloaded as a zip file
from GitHub and is installed under 'D:\\USERDATA\\pylib\\', your folder should now
resemble the screenshot below.

.. image:: pylib.png

.. _here: https://github.com/ajwdewit/pcse/archive/master.zip


Testing the PCSE package
========================

To guarantee its integrity, the PCSE package includes a number of self
tests that test individual components as well as the entire simulation. These tests
verify that the output produced by the different components matches with the
expected outputs. Test data for the individual components can be found
in the `pcse.tests.test_data` package, while the test data for the entire chain
is stored in an SQLite database (pcse.db). This database can be found under
`.pcse` in your home folder and will be automatically created when importing
PCSE for the first time. When you delete the database file manually it will be
recreated next time you import PCSE.

For importing PCSE we first need to add the location of PCSE ('D:\\USERDATA\\pylib\\')
to the search path of Python::

    C:\>python
    Python 2.7.9 |Continuum Analytics, Inc.| (default, Dec 18 2014, 17:00:07) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://binstar.org    Enthought Python Distribution -- www.enthought.com
    >>> import sys
    >>> sys.path.append(r"D:\USERDATA\pylib")

Next, PCSE can be imported and the tests can be executed by calling
the `test()` function at the top of the package::

    >>> import pcse
    Building PCSE demo database at: C:\Users\wit015\.pcse\pcse.db
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
due to problems with the test. Moreover, SQLAlchemy may complain with a
warning that can be safely ignored::

     /usr/lib/python2.7/dist-packages/sqlalchemy/types.py:307: SAWarning:
     Dialect sqlite+pysqlite does *not* support Decimal objects natively, and
     SQLAlchemy must convert from floating point - rounding errors and other
     issues may occur. Please consider storing Decimal numbers as strings or
     integers on this platform for lossless storage.
         d[coltype] = rp = d['impl'].result_processor(dialect, coltype)

