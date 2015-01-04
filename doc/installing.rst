***************
Installing PCSE
***************

Requirements and dependencies
=============================

PCSE is being developed on Ubuntu Linux 12.04 using python 2.7 and is known to work with
the 3.x series (using the 2to3 tool). As python is a platform independent language, PCSE
works equally well on Windows or Mac OSX.  The most straightforward approach for installing
python is through one of the prepackaged python distributions such as `Enthought Canopy`_,
`Anaconda`_ or `PythonXY`_. The following screen dump shows the version of python, numpy and
SQLAlchemy that were used to develop PCSE::

    Python 2.7.6 (default, Dec 16 2013, 12:39:22)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy as np
    >>> np.__version__
    '1.8.0'
    >>> import sqlalchemy as sa
    >>> sa.__version__
    '0.8.4'

.. _Enthought Canopy: https://www.enthought.com/products/canopy/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _PythonXY: https://code.google.com/p/pythonxy/wiki/Welcome

How to install
==============

PCSE can be installed in different ways and the best way to install it depends on your
requirements. The most convenient option to install PCSE is through the Python Package
Index (PYPI). Installing from PYPI is mostly useful if you are interested in using the functionality
provided by PCSE in your own scripts, but are not interested in modifying or contributing to
PCSE itself. Installing from PyPI is done using the package installer `pip`, but this
will only work when you have write access on the python site-packages
folder. If you do not have this permission you should install PCSE into a
`virtual environment`_.::

    allard$ pip install PCSE
    Downloading/unpacking PCSE
      Running setup.py install for PCSE
      ...
      ...
    Successfully installed PCSE
    Cleaning up...
    allard$

.. _virtual environment: http://docs.python-guide.org/en/latest/dev/virtualenvs/

The second option is useful if you want to develop or contribute to PCSE.
In that case you should fork my `PCSE
repository`_ on github and get a local copy of PCSE using `git clone`. See the help on github_
and for windows user the `github for windows`_

.. _github for windows: https://windows.github.com/
.. _github: https://help.github.com/
.. _PCSE repository: https://github.com/ajwdewit/pcse

Finally, the last option is to download the PCSE package as a zip file from GitHub.com
using the link `here`_. Just unzip the package at a suitable location.
Note that the top directory in the zip file is `pcse-<branchname>`.
The actual PCSE is package is inside this folder and needs to be put on your file system.
For running the test included with PCSE, we assume that PCSE was downloaded as a zip file
from github and is installed under 'D:\\USERDATA\\pylib\\', your folder should now
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
to the search path of python::

    C:\>python
    Enthought Python Distribution -- www.enthought.com
    Version: 7.0-2 (32-bit)

    Python 2.7.1 |EPD 7.0-2 (32-bit)| (r271:86832, Dec  2 2010, 10:35:02) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
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
    runTest (pcse.tests.test_wofost.TestWaterlimitedPotato) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterWheat) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedGrainMaize) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedWinterRapeseed) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialPotato) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSpringBarley) ... ok
    runTest (pcse.tests.test_wofost.TestWaterlimitedSunflower) ... ok
    runTest (pcse.tests.test_wofost.TestPotentialWinterWheat) ... ok

    ----------------------------------------------------------------------
    Ran 19 tests in 29.748s

    OK
    >>>

If the model output matches the expected output the test will report 'OK',
otherwise an error will be produced with a detailed traceback on where the
problem occurred. Note that the results may deviate from the output above
because one or more
tests may have been temporarily disabled (skipped) often due to problems
with the test. Moreover, SQLAlchemy may complain with a warning that can be safely ignored::

     /usr/lib/python2.7/dist-packages/sqlalchemy/types.py:307: SAWarning:
     Dialect sqlite+pysqlite does *not* support Decimal objects natively, and
     SQLAlchemy must convert from floating point - rounding errors and other
     issues may occur. Please consider storing Decimal numbers as strings or
     integers on this platform for lossless storage.
         d[coltype] = rp = d['impl'].result_processor(dialect, coltype)

