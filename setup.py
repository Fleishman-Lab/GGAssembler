#!/usr/bin/env python

# from distutils.extension import Extension

from distutils.command.build import build as build_orig

import pkg_resources
from setuptools import Extension, find_packages, setup


class build(build_orig):
    def finalize_options(self):
        super().finalize_options()
        try:
            __builtins__.__NUMPY_SETUP__ = False
        except AttributeError:
            pass
        ext = next(m for m in self.distribution.ext_modules if m == ext_modules[0])
        ext.include_dirs.append(
            pkg_resources.resource_filename("numpy", "core/include")
        )


#
# directive_defaults = Options.get_directive_defaults()
# directive_defaults['profile'] = True
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True

USE_CYTHON = "auto"

if USE_CYTHON:
    try:
        from Cython.Build import cythonize
        from Cython.Compiler import Options

        Options.docstrings = True
    except ImportError:
        if USE_CYTHON == "auto":
            USE_CYTHON = False
        else:
            raise

ext = ".pyx" if USE_CYTHON else ".cpp"

ext_modules = [
    Extension(
        "dawdlib.dijkstra.colorful",
        ["dawdlib/dijkstra/colorful" + ext],
        include_dirs=[],
        language="c++",
        extra_compile_args=["-Ofast", "-std=c++11"]
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
    ext_modules=ext_modules,
    cmdclass={"build": build},
    install_requires=[
        "networkx",
        "numpy",
        "biopython",
        "synbiochem-py",
        "pandas",
        "fire",
        "tabulate",
    ],
    setup_requires=["numpy"],
    include_package_data=True,
    entry_points={"console_scripts": ["dawdlib=dawdlib.__main__:main"]},
)
