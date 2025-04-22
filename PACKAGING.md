# GenomeKit Packaging Guide

This document explains how to build, install, and distribute the GenomeKit package.

## Package Structure

The package has been organized according to standard Python package conventions:

```
genomekit/                 # Main package directory
├── __init__.py            # Package initialization
├── cli.py                 # Command-line interface
└── modules/               # Module directory
    ├── __init__.py        # Module initialization
    ├── microarray_generator.py
    ├── opencravat_analyzer.py
    └── snp_checker.py
```

## Building the Package

To build the package, use these commands:

```bash
# Install build tools if not already installed
pip install build wheel

# Build the package
python -m build
```

This will create both source distribution and wheel packages in the `dist/` directory.

## Installing the Package

### Development Installation

For development, you can install the package in "editable" mode:

```bash
# Install in development mode
pip install -e .
```

This allows you to modify the code and see changes immediately without reinstalling.

### Regular Installation

To install from the built packages:

```bash
pip install dist/genomekit-0.1.0-py3-none-any.whl
```

Or from source directory:

```bash
pip install .
```

## Using the Package

After installation, you can use the command-line interface:

```bash
# Show available commands
genomekit --help

# Run SNP checker
genomekit snps path/to/file.vcf --section metabolism

# Generate microarray files
genomekit microarray --bam path/to/file.bam --formats 23andMe_V5 --outdir ./output
```

Or import the modules in your Python scripts:

```python
from genomekit.modules.snp_checker import SNPChecker
from genomekit.modules.microarray_generator import MicroarrayGenerator

# See examples/ directory for more detailed examples
```

## Distributing the Package

If you want to distribute the package publicly:

1. Register an account on PyPI (https://pypi.org/)
2. Install twine: `pip install twine`
3. Upload your package: `twine upload dist/*`

For private distribution, you can share the wheel file directly or set up a private PyPI server.

## Dependencies

The package has the following dependencies:
- pandas
- numpy
- biopython
- pysam

External tools that should be installed on the system:
- bcftools
- tabix
- GATK (optional)
