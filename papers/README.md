# Reproducing Numerical Results

## Overview
This directory provides the scripts and instructions required to reproduce the numerical results reported in the paper.
The workflow is based on the **PySCF** package, which must be installed before running any of the scripts.

## Quick Start
Install **PySCF** using either `pip` or `conda`:

```bash
# Using pip
pip install pyscf

# Or, using conda
conda install -c conda-forge pyscf
```

## Data Preparation

Before performing calculations, you need to generate the FCIDUMP files:

1. Navigate to the `../data` directory.
1. Run:
```bash
python3 generate_fcidump.py
```
This will produce several files with the `.fcidump` extension.

> **Note**: These `.fcidump` files can be large and are therefore not included in the Git repository.
> Make sure you have sufficient disk space available before generating them.

## Running the Calculations

Once the `.fcidump` files are ready, you can execute the scripts in this directory to reproduce the results presented in the paper.
The script `run_all.sh` contains one line per experiment, each corresponding to a specific calculation.