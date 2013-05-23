#!/usr/bin/env python

# NOTE: multiprocessing import required for issues with nose tests.
# See: http://bugs.python.org/issue15881#msg170215
import multiprocessing

from setuptools import setup

setup(
    name='optmage',
    version='0.9',
    author='Church Lab',
    author_email='gleb@mit.edu',
    maintainer='Gleb Kuznetsov',
    maintainer_email='gleb@mit.edu',
    url='http://arep.med.harvard.edu/optMAGE/',
    package_dir={'': 'src'},
    packages = ['optmage'],
    install_requires=['biopython >= 1.6.1'],
    test_suite = 'nose.collector',
)
