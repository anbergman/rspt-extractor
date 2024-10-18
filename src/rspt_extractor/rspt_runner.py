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

import argparse
import numpy as np
from rspt_extractor import (
    RsptScf,
    RsptExchange,
    downscale_exchange,
    extract_projections,
    print_projections,
)


# def default_workflow(xstr="100", ystr="010", zstr="001"):
def default_workflow(input_directions=None, maptype="C", atomlist=None):
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
    # Define input directions and atoms
    if input_directions is None:
        input_directions = ["100", "010", "001"]
    input_atoms = range(1, 5)

    exchange_data: dict[str, list[RsptExchange]] = {}

    # Extract information from each input file
    for direction in input_directions:
        print("Extracting exchange data from:", direction)
        exchange_data[direction] = []
        for atoms in input_atoms:
            # print(f"Atom: {atoms}")
            file_name = f"spin-{direction}/out-{atoms}"
            extracted_data = RsptExchange(file_name, maptype)
            exchange_data[direction].append(extracted_data.outmap)

    concatenated_data: dict[str, np.ndarray] = {}
    for key, value in exchange_data.items():
        concatenated_data[key] = np.concatenate([v for v in value])

    # Masks for extracting the "best" estimates for J, D, A
    # Example: spin-axis x gives Dx and Ax but J as 0.5*(Jyy*Jzz)
    # This should not be system dependent.
    x_mask = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 1.0], [0.0, 1.0, 0.5]])
    y_mask = np.array([[0.5, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.5]])
    z_mask = np.array([[0.5, 1.0, 0.0], [1.0, 0.5, 0.0], [0.0, 0.0, 0.0]])
    mask_list = [x_mask, y_mask, z_mask]

    masked_data: dict[str, np.ndarray] = {
        key: value.copy() for key, value in concatenated_data.items()
    }

    for idx, (key, value) in enumerate(masked_data.items()):
        for row in value:
            j_mat = row[5:14].reshape(3, 3)
            j_mat = mask_list[idx] * j_mat * 3.0
            row[5:14] = j_mat.flatten()

    summed_data = (
        np.array(masked_data["001"] + masked_data["010"] + masked_data["100"]) / 3.0
    )

    # Filter out the Jij values so that only selected atoms are included
    # eg. only the magnetic atoms
    # j_truncated = downscale_exchange(summed_data, [1, 2])
    print("Atomlist:", atomlist)
    if atomlist is None:
        j_truncated = np.array(summed_data, dtype=np.float32)
    else:
        j_truncated = downscale_exchange(summed_data, atomlist)

    # Create the various projections you can make from the exchange interactions
    j_dict = extract_projections(j_truncated)
    # Print all projections to file
    print_projections(j_dict, maptype)

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


def run_scf(filepath):
    """
    Run the SCF (Self-Consistent Field) calculation using the given file path.

    This function initializes an RsptScf object with the provided file path and
    generates several output files related to the lattice, atomic positions,
    magnetic moments, and a template for further calculations.

    Args:
        filepath (str): The path to the input file required for the SCF
                        calculation.

    Outputs:
        lattice.dat: File containing lattice information.
        posfile: File containing atomic positions.
        momfile: File containing magnetic moments.
        inpsd.minimal: Template file for further calculations.

    Author:
        Anders Bergman
    """
    rspt_exchange = RsptScf(filepath)
    rspt_exchange.print_lattice("lattice.dat")
    rspt_exchange.print_positions("posfile")
    rspt_exchange.print_moments("momfile")
    rspt_exchange.print_template("posfile", "momfile", "j_scalar.dat", "inpsd.minimal")


def run_exchange(filepath, maptype="C"):
    """
    Executes the RSPT exchange process and saves the output to a file.

    This function initializes an RsptExchange object with the provided
    filepath, runs the exchange process, and saves the output to a file
    named 'jfile.tensor'.

    Args:
        filepath (str): The path to the input file for the RSPT exchange
        process.

    Returns:
        None

    Author:
        Anders Bergman
    """
    rspt_exchange = RsptExchange(filepath, maptype)
    rspt_exchange.save_output("jfile.tensor")


def main():
    """
    Main function to parse command-line arguments and run the appropriate
    workflow based on the provided options.

    This function sets up a command-line interface (CLI) with three options:
    - `-s` or `--scf`: Run the SCF workflow with the given file path.
    - `-e` or `--exchange`: Run the exchange workflow with the given file path.
    - `-r` or `--run`: Run the default workflow with the given directories.

    The function then parses the arguments and calls the corresponding
    workflow function based on the provided options.

    Args:
        None

    Returns:
        None

    Raises:
        SystemExit: If the argument parsing fails or if an invalid combination
        of arguments is provided.

    Example:
        To run the SCF workflow:
            $ python rspt_runner.py -s /path/to/scf_file

        To run the exchange workflow:
            $ python rspt_runner.py -e /path/to/exchange_file

        To run the default workflow:
            $ python rspt_runner.py -r dirx diry dirz

    Author:
        Anders Bergman
    """
    parser = argparse.ArgumentParser(description="RSPt exchange extraction tool")

    # CLI options
    parser.add_argument(
        "-s", "--scf", type=str, help="Run SCF workflow with the given FILE_PATH"
    )
    parser.add_argument(
        "-e",
        "--exchange",
        type=str,
        help="Run exchange workflow with the given FILE_PATH",
    )
    parser.add_argument(
        "-r",
        "--run",
        nargs=3,
        metavar=("DIRX", "DIRY", "DIRZ"),
        help="Run default workflow with directories",
    )
    parser.add_argument(
        "-m",
        "--maptype",
        type=str,
        help="Output maptype: (C)artesian, (D)irect or maptype (3)",
    )
    parser.add_argument(
        "-a",
        "--atoms",
        nargs="+",
        type=int,
        help="List of atoms to extract exchange interactions for",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Decide which function to run based on the CLI arguments
    if args.scf:
        run_scf(args.scf)
    elif args.exchange:
        run_exchange(args.exchange, args.maptype)
    elif args.run:
        default_workflow(
            input_directions=args.run, maptype=args.maptype, atomlist=args.atoms
        )
        # default_workflow(args.run[0], args.run[1], args.run[2])
    else:
        default_workflow(
            input_directions=None, maptype=args.maptype, atomlist=args.atoms
        )


if __name__ == "__main__":
    main()
