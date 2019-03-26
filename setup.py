from __future__ import print_function
from setuptools import setup, find_packages
import os
import io

PACKAGE = "pcse"
NAME = "PCSE"
DESCRIPTION = 'Framework for developing crop simulation models, ' \
              'includes an implementation of '\
              'the WOFOST and LINTUL crop simulation models.'
AUTHOR = "Allard de Wit"
AUTHOR_EMAIL = 'allard.dewit@wur.nl'
URL = 'http://github.com/ajwdewit/pcse/'
LICENSE="EUPL"
VERSION = "5.4.2"

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


long_description = read('README.rst')

setup(
    name=NAME,
    version=VERSION,
    url=URL,
    download_url='http://github.com/ajwdewit/pcse/tarball/'+VERSION,
    license='EUPL',
    author=AUTHOR,
    install_requires=['numpy>=1.6.0',
                      'SQLAlchemy>=0.8.0',
                      'PyYAML>=3.11',
                      'xlrd>=0.9.3',
                      'xlwt>=1.0.0',
                      'requests>=2.0.0',
                      'pandas>=0.20',
                      'traitlets-pcse==5.0.0.dev'],
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    platforms='any',
    test_suite='pcse.tests.make_test_suite',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.1 (EUPL 1.1)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering']
)

