
# -- import packages: ----------------------------------------------------------
import setuptools
import re
import os
import sys


# -- run setup: ----------------------------------------------------------------
setuptools.setup(
    name="adata_query",
    version="0.0.1rc1",
    python_requires=">3.9.0",
    author="Michael E. Vinyard",
    author_email="mvinyard.ai@gmail.com",
    url="https://github.com/mvinyard/AnnDataQuery",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="Fetch and format data matrices from AnnData.",
    packages=setuptools.find_packages(),
    install_requires=[
        "anndata>=0.9.1",
        "torch>=2.0.1",
        "autodevice>=0.0.2",
        "ABCParse>=0.0.6",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
