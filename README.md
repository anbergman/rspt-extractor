# RSPt Extractor

**Version:** 0.1.0
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
