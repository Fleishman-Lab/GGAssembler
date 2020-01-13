#!/usr/bin/env python

from setuptools import setup, find_packages
from distutils.extension import Extension
import numpy

# import Cython
from Cython.Compiler import Options
Options.docstrings = True
#
# directive_defaults = Options.get_directive_defaults()
# directive_defaults['profile'] = True
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True

USE_CYTHON = True

if USE_CYTHON:
    try:
        from Cython.Build import cythonize
    except ImportError:
        if USE_CYTHON == 'auto':
            USE_CYTHON = False
        else:
            raise

ext = '.pyx' if USE_CYTHON else '.c'

ext_modules = [
    Extension(
        'dawdlib.dijkstra.colorful',
        ['dawdlib/dijkstra/colorful'+ext],
        include_dirs=[numpy.get_include()],
        # define_macros=[('CYTHON_TRACE', '1'), ('CYTHON_TRACE_NOGIL', '1')]
    )
]

if USE_CYTHON:
    ext_modules = cythonize(ext_modules, gdb_debug=False)


setup(
    name="dawdlib",
    version="0.0.1a0",
    url="https://github.com/Fleishman-Lab/dawdlib.git",
    description="Utils to help create a cost sensitive degenerate codon sequence",
    packages=find_packages(),
    zip_safe=False,
    ext_modules=ext_modules
)
