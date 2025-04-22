#!/usr/bin/env python3
# coding: utf-8

from setuptools import setup, find_packages

# Read long description from README.md
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="genome",
    version="0.1.0",
    author="Alex Skryl",
    author_email="rut216@gmail.com",
    description="A toolkit for genome analysis and microarray generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/skryl/genome",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "numpy",
        "biopython",
        "pysam",
    ],
    entry_points={
        "console_scripts": [
            "genomekit=genomekit.cli:main",
        ],
    },
)
