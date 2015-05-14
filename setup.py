from __future__ import print_function
from setuptools import setup, find_packages
import os
import io

#import pcse

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
    name='PCSE',
    #version=pcse.__version__,
    version="5.1.1",
    url='http://github.com/ajwdewit/pcse/',
    download_url='http://github.com/ajwdewit/pcse/tarball/5.1.1',
    license='EUPL',
    author='Allard de Wit',
    install_requires=['numpy>=1.6.0',
                      'SQLAlchemy>=0.8.0',
                      'xlrd>0.9.0',
                      'tabulate>=0.7.0'],
    author_email='allard.dewit@wur.nl',
    description='Framework for developing crop simulation models, includes an '
                'implementation of the WOFOST crop simulation model.',
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
#    package_data = {
#        # Include any files *.txt, *.rst, *.conf, *.csv or *.sql files:
#        '': ['*.txt', '*.rst', '*.conf', '*.csv', '*.sql'],
#    },
    platforms='any',
    test_suite='pcse.tests.make_test_suite',
    use_2to3=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.1 (EUPL 1.1)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Fortran',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering']
)

