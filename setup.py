# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

from setuptools import setup, find_packages
from setuptools.dist import Distribution
from pathlib import Path

#######################################
# Prepare list of compiled extensions #
#######################################

extensions = []


class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name."""
    def has_ext_modules(_):  # noqa
        return True


#########
# Setup #
#########

version_file = Path(__file__).parent / 'xtrack/_version.py'
dd = {}
with open(version_file.absolute(), 'r') as fp:
    exec(fp.read(), dd)
__version__ = dd['__version__']

setup(
    name='xtrack',
    version=__version__,
    description='Tracking library for particle accelerators',
    long_description='Tracking library for particle accelerators',
    url='https://xsuite.readthedocs.io/',
    author='G. Iadarola et al.',
    license='Apache 2.0',
    download_url="https://pypi.python.org/pypi/xtrack",
    project_urls={
        "Bug Tracker": "https://github.com/xsuite/xsuite/issues",
        "Documentation": 'https://xsuite.readthedocs.io/',
        "Source Code": "https://github.com/xsuite/xtrack",
    },
    packages=find_packages(),
    ext_modules=extensions,
    install_requires=[
        'numpy>=1.0',
        "pandas>=2.0",
        'scipy',
        'tqdm',
        'xobjects',
        'xpart',
        'xdeps'
    ],
    extras_require={
        'tests': ['cpymad', 'nafflib', 'PyHEADTAIL', 'pytest', 'pytest-mock'],
    },
    # The very non-obvious way files can be included in either sdist or bdist.
    # List the file in:
    # - MANIFEST.in if you want it in both sdist and bdist,
    # - package_data and exclude from the manifest if you want it in bdist only,
    # - exclude_package_data if you want it in sdist only,
    # - nowhere if you don't want to package them.
    # The following enables us to publish a source-only distribution and binary
    # wheels with the prebuilt kernels in the following way:
    # 1. Source
    #   1.1. Build sdist as usual
    # 2. Binary
    #   2.1. Install the package with '-e'.
    #   2.2. Prebuild the kernels, make sure they are in xtrack/prebuilt_kernels.
    #   2.3. Build the bdist wheel.
    include_package_data=True,
    package_data={
        'xtrack': [
            'prebuilt_kernels/*.so',
            'prebuilt_kernels/*.dylib',
            'prebuilt_kernels/*.dll',
            'prebuilt_kernels/*.json',
        ]
    },
    distclass=BinaryDistribution,
)
