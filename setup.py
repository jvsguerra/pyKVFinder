# System imports
# from distutils.core import setup, Extension
from setuptools import setup, Extension, dist

# pyKVFinder information
from pyKVFinder import _name, _version
# _name = "pyKVFinder"
# _version = "0.1"

# Prepare reqs from requirements.txt
with open('requirements.txt') as f:
    reqs = f.read().splitlines()

# Third-party modules - we depend on numpy for everything
dist.Distribution().fetch_build_eggs([req for req in reqs if req.find('numpy') != -1])
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Extension modules
_grid = Extension(
    name="_gridprocessing",
    sources=["src/grid.i", "src/grid.c"],
    include_dirs=[numpy_include],
    extra_compile_args=['-fopenmp', '-lm'],
    extra_link_args=['-lgomp'],
    # swig_opts=['-c']
)

# Setup
setup(
    name=_name, 
    version=_version,
    ext_modules=[_grid],
    include_package_data=True,
    install_requires=reqs,
    packages=['pyKVFinder'],
    entry_points={'console_scripts': ['pyKVFinder=pyKVFinder.main:run']},
    # extra_compile_args=['-fopenmp'],
    # extra_link_args=['-lgomp'],
    # swig_opts=['-threads']
)
