"""
Module for extracting and processing exchange data from RSPt calculations.

This module reads exchange data from multiple input files, processes the data
by applying masks, and then combines the data. It also extracts and prints
information from an SCF run.

Author: Anders Bergman

Functions:
    main(): Main function to execute the data extraction and processing.

Usage:
    This script is intended to be run as a standalone module.

Example:
    $ python rspt_runner.py

Dependencies:
    - numpy
    - rspt_scf
    - rspt_exchange (as rspt_xc)

Attributes:
    input_directions (list): List of input directions.
    input_atoms (range): Range of input atoms.
    input_files (list): List to store input files.
    exchange_data (dict): Dictionary to store extracted exchange data.
    concatenated_data (dict): Dictionary to store concatenated exchange data.
    x_mask (np.array): Mask for x direction.
    y_mask (np.array): Mask for y direction.
    z_mask (np.array): Mask for z direction.
    mask_list (list): List of masks for each direction.
    masked_data (dict): Dictionary to store masked exchange data.
    summed_data (np.array): Array to store summed exchange data.
    j_truncated (np.array): Array to store downscaled exchange data.
    j_dict (dict): Dictionary to store extracted projections.
    scf_file (str): File path for SCF data.
    scf_data (RsptScf): Object to store SCF data.

"""

import numpy as np
from rspt_extractor import (
    RsptScf,
    RsptExchange,
    downscale_exchange,
    extract_projections,
    print_projections,
)


def main():
    """
    Main function to extract and process exchange data from RSPT output files.

    This function performs the following steps:
    1. Defines input directions and atoms.
    2. Extracts exchange data from specified input files.
    3. Concatenates the extracted data.
    4. Applies masks to the concatenated data.
    5. Computes the average of the masked data.
    6. Downscales the exchange data.
    7. Extracts and prints projections from the downscaled data.
    8. Extracts and prints SCF data from the last input direction.
    9. Prints lattice, basis, and moments information.
    10. Saves lattice, positions, moments, and template data to files.

    Args:
        None

    Returns:
        None

    Raises:
        FileNotFoundError: If any of the input files are not found.
        ValueError: If there is an issue with data extraction or processing.

    Author:
        Anders Bergman
    """
    # List of input files
    input_directions = ["100", "010", "001"]
    input_atoms = range(1, 5)

    # input_files: list[str] = []

    exchange_data: dict[str, list[RsptExchange]] = {}

    # Extract information from each input file
    for direction in input_directions:
        print("Extracting exchange data from:", direction)
        exchange_data[direction] = []
        for atoms in input_atoms:
            print(f"Atom: {atoms}")
            file_name = f"spin-{direction}/out-{atoms}"
            extracted_data = RsptExchange(file_name)
            exchange_data[direction].append(extracted_data.outmap)

    concatenated_data: dict[str, np.ndarray] = {}
    # concatenated_data = {}
    for key, value in exchange_data.items():
        print(key, len(value))
        concatenated_data[key] = np.concatenate([v for v in value])

    x_mask = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 1.0], [0.0, 1.0, 0.5]])
    y_mask = np.array([[0.5, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.5]])
    z_mask = np.array([[0.5, 1.0, 0.0], [1.0, 0.5, 0.0], [0.0, 0.0, 0.0]])
    mask_list = [x_mask, y_mask, z_mask]

    # masked_data = {key: value.copy() for key, value in concatenated_data.items()}
    masked_data: dict[str, np.ndarray] = {
        key: value.copy() for key, value in concatenated_data.items()
    }
    # masked_data = {
    #     key: [rspt_xc.RsptExchange(v) for v in value]
    #     for key, value in exchange_data.items()
    # }

    for idx, (key, value) in enumerate(masked_data.items()):
        for row in value:
            j_mat = row[5:14].reshape(3, 3)
            j_mat = mask_list[idx] * j_mat * 3.0
            row[5:14] = j_mat.flatten()

    summed_data = (
        np.array(masked_data["001"] + masked_data["010"] + masked_data["100"]) / 3.0
    )

    j_truncated = downscale_exchange(summed_data, [1, 2, 3, 4])
    j_dict = extract_projections(j_truncated)
    print_projections(j_dict)

    # Exctract data from an scf-run
    scf_file = f"spin-{input_directions[-1]}/out-scf"
    print("Extracting SCF data from:", scf_file)
    scf_data = RsptScf(scf_file)

    # Print the lattice file from any of the input files
    print("Lattice:")
    print(scf_data.lattice)
    print("Basis:")
    print(scf_data.basis)
    print("Moments:")
    print(scf_data.moments)
    scf_data.print_lattice("lattice.dat")
    scf_data.print_positions("posfile")
    scf_data.print_moments("momfile")
    scf_data.print_template("posfile", "momfile", "j_scalar.dat", "inpsd.minimal")
    print("Extraction and storage completed successfully!")


if __name__ == "__main__":
    main()
