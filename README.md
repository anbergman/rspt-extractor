# RSPt Extractor

**Version:** 0.3.0
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

The classes included in the package can be customized to your preferred directory structure with the suitable flags. Files needed from an `RSPt` calculation are the following:

- `data`: `RSPt` structure input file. Needed for positions of atoms.
- `out.scf` (arbitrary name): Terminal output from `RSPt` SCF calculation. Needed for magnetic moments.
- `out.jij` (arbitrary name): Terminal output from `RSPt` Jij calculation. Needed for extraction of exchange interactions and neighbour information.

To allow for flexible naming and placing of these files, all inputs are given as flags to the parser. Thus, assuming the naming convention above, and a single Jij output file, the syntax would be

```bash
rspt-parser --data data --scf out.scf --exchange out.jij
```

to extract the calculated exchange interactions.

## Flags

While the package can be customized to your own needs, the default script `rspt-extractor` comes with a number of flags to allow more user control. The flags are visible by `rspt-extractor --help` and are currently as follows

```syntax
RSPt exchange extraction tool

options:
  -h, --help            show this help message and exit
  -d DATA, --data DATA  FILE_NAME for RSPt input data
  -s SCF, --scf SCF     FILE_NAME for RSPt SCF output data
  -e FILE_NAME [FILE_NAME ...], --exchange FILE_NAME [FILE_NAME ...]
                        Run exchange workflow for the given FILE_NAMEs
  --rx FILE_NAME [FILE_NAME ...]
                        Input files for relativistic exchange workflow (x-direction)
  --ry FILE_NAME [FILE_NAME ...]
                        Input files for relativistic exchange workflow (y-direction)
  --rz FILE_NAME [FILE_NAME ...]
                        Input files for relativistic exchange workflow (z-direction)
  -m MAPTYPE, --maptype MAPTYPE
                        Output maptype: (C)artesian, (D)irect or maptype (3)
  -a ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]
                        List of atoms to extract exchange interactions for
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for moment magnitudes. Default is 0.0
  -c CUTOFF, --cutoff CUTOFF
                        Radius for excange interaction cutoff
```

## Included examples

### FePt

An example for relativistic exchange interaction extraction can be found in the folder `examples/FePt`. For relativistic interactions, *three* different calculations with different spin axes are needed. In this example the calculations are grouped within the following folder structure

```bash
├── FePt
│   ├── spin-001
│   ├── spin-010
│   └── spin-100
```

with the output files from both SCF and Jij calculations available in each subfolder. To extract the relativistic exchange interactions the syntax for this example should be

```syntax
rspt-parser --scf spin-001/out-scf --data spin-001/data --rx spin-100/out-? --ry spin-010/out-? --rz spin-001/out-?
```

where `--rx`, `--ry`, and `--rz` controls where the Jij output files are fund. Notice the wild-carding since there are several Jij output files in each directory.

*NOTE:* The relativistic workflow is not guaranteed to work due to the risk of uncontrolled reshuffeling of atoms between the different spin-axis calculations. A fix for this is pending for the `RSPt` code.

## Additional flags

Thresholding the output can be done by either moment magnitude or cutoff radius. To only extract the Fe atoms in the `FePt` example above we can add the `--threshold 0.5` which will filter out all atoms with a magnetic moment lower than $0.5 \mu_B$. To enforce a cutoff radius, the flag `--cutoff 3.0` can be included to truncate the exchange interactions to only include interactions within $3.0$ lattice parameters. Example:

```syntax
rspt-parser --scf spin-001/out-scf --data spin-001/data --rx spin-100/out-? --ry spin-010/out-? --rz spin-001/out-? --threshold 0.5 --cutoff 3.0
```
