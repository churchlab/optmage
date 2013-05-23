#!/usr/bin/env python

# NOTE: multiprocessing import required for issues with nose tests.
# See: http://bugs.python.org/issue15881#msg170215
import multiprocessing

from setuptools import setup

setup(
    name='foo',
    version='0.9',
    py_modules=['foo'],
    author='Church Lab',
    author_email='gleb@mit.edu',
    maintainer='Gleb Kuznetsov',
    maintainer_email='gleb@mit.edu',
    url='http://arep.med.harvard.edu/optMAGE/',
    package_dir={'': 'src'},
    packages = ['optmage'],
    test_suite = 'nose.collector',
)
