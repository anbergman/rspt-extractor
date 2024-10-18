# RSPt Extractor

**Version:** 0.2.0
**Author:** Anders Bergman

## Overview

`rspt-extractor` is a Python package for extracting and processing data from **RSPt calculations**. It reads exchange data from multiple input files, processes the data, and combines the results. The package also provides methods to extract and print information from an **SCF** (Self-consistent field) run.

This tool is specifically designed to help researchers working with **RSPt** and related calculations in condensed matter physics.

## Features

- **Exchange Data Extraction**: Extracts exchange data from multiple files and processes it.
- **Data Masking**: Applies masks to filter the data based on user-specified criteria.
- **SCF Run Information**: Extracts and prints key details from an SCF calculation.

## Installation

To install `rspt-extractor`, use the following command after cloning the repository:

```bash
pip install -e .
```

## Usage

The classes included in the package can be customized to your preferred directory structure. Assuming the default structure for relativistic calculations i.e. three folders named `spin-001`, `spin-010`, and `spin-100`, then the built-in script `rspt-runner` can be invoked as

```bash
rspt-runner
```

This should work for the included examples.

## Options

While the package can be customized to your own needs, the default script `rspt-runner` comes with a number of flags to allow more user control. The flags are visible by `rspt-runner --help` and are currently as follows

```syntax
RSPt exchange extraction tool

options:
  -h, --help            show this help message and exit
  -s SCF_PATH, --scf SCF_PATH     Run SCF workflow with the given SCF_PATH
  -e XC_PATH, --exchange XC_PATH
                        Run exchange workflow with the given XC_PATH
  -r DIRX DIRY DIRZ, --run DIRX DIRY DIRZ
                        Run default workflow with directories DIRX DIRY DIRZ
  -m MAPTYPE, --maptype MAPTYPE
                        Output maptype: (C)artesian, (D)irect or maptype (3)
  -a ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]
                        List of atoms to extract exchange interactions for. Default all.
```
