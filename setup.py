#!/usr/bin/env python

"""Setup script for ssw-aligner.

Extracted StripedSmithWaterman (SSW) C implementation and Python/Cython wrapper.

Original SSW C library:
  Copyright (c) 2012-2015 Boston College. MIT License.

Original Cython wrapper:
  Copyright (c) 2013--, scikit-bio development team. Modified BSD License.
"""

import platform
import os
import subprocess
import sysconfig

from setuptools import setup, find_packages
from setuptools.extension import Extension

import numpy as np
from Cython.Build import cythonize


# --- Compiler detection (adapted from scikit-bio) ---

def check_bin(ccbin, source, allow_dash):
    """Check if a given compiler matches the specified name."""
    source0 = source.split()[0]
    bsource = os.path.basename(source0)
    if allow_dash:
        found = False
        for el in bsource.split("-"):
            if el == ccbin:
                found = True
                break
    else:
        found = bsource == ccbin
    return found


clang = False
icc = False
gcc = True

try:
    if os.environ["CC"] == "gcc":
        gcc = True
    elif os.environ["CC"] != "":
        gcc = False
except KeyError:
    pass

if not gcc:
    try:
        if check_bin("clang", os.environ["CC"], False):
            clang = True
        elif check_bin("icc", os.environ["CC"], True):
            icc = True
    except KeyError:
        pass
else:
    try:
        if check_bin("clang", sysconfig.get_config_vars()["CC"], False):
            clang = True
            gcc = False
        elif check_bin("icc", sysconfig.get_config_vars()["CC"], True):
            icc = True
            gcc = False
    except KeyError:
        pass

if gcc:
    try:
        if (
            subprocess.check_output(["gcc", "--version"], universal_newlines=True).find("clang") != -1
        ):
            clang = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

# --- Compile flags ---

ssw_extra_compile_args = ["-I."]

if platform.system() != "Windows":
    if icc:
        ssw_extra_compile_args.extend(["-qopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
    elif not clang:
        ssw_extra_compile_args.extend(["-fopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
elif platform.system() == "Windows":
    ssw_extra_compile_args.extend(["-openmp:experimental"])

if platform.machine() == "i686":
    ssw_extra_compile_args.append("-msse2")

# --- Extensions ---

ext = ".pyx"
extensions = [
    Extension(
        "ssw_aligner._ssw_wrapper",
        [
            "ssw_aligner/_ssw_wrapper" + ext,
            "ssw_aligner/_lib/ssw.c",
        ],
        extra_compile_args=ssw_extra_compile_args,
        include_dirs=[np.get_include()],
    ),
]

extensions = cythonize(extensions, force=True)

# --- Setup ---

setup(
    name="ssw-aligner",
    version="0.1.0",
    license="BSD-3-Clause AND MIT",
    description="Standalone StripedSmithWaterman (SSW) aligner - extracted from scikit-bio.",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="scikit-bio development team / Boston College (SSW C library)",
    url="https://github.com/scikit-bio/scikit-bio",
    packages=find_packages(),
    package_data={"ssw_aligner": ["py.typed", "_ssw_wrapper.pyi"]},
    ext_modules=extensions,
    include_dirs=[np.get_include()],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.17.0",
    ],
    extras_require={
        "dev": ["pytest", "cython", "scikit-bio==0.6.2"],
    },
)
