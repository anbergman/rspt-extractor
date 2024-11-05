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
    scf_data (RsptData): Object to store SCF data.

"""

import argparse
import numpy as np
from rspt_extractor import (
    RsptData,
    RsptExchange,
    downscale_exchange,
    extract_projections,
    print_projections,
    extract_position_data,
)


def rspt_get_scf_data(args):
    """
    Read RSPt SCF (Self-Consistent Field) data using the given file path.

    This function initializes an RsptData object with the provided file path and
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
    rspt_data = RsptData(args)
    if rspt_data.scf_file:
        rspt_data.print_lattice("lattice.dat")
    if rspt_data.data_file:
        rspt_data.print_positions("posfile")
    if rspt_data.data_file and rspt_data.scf_file:
        rspt_data.print_moments("momfile")
    rspt_data.print_template("posfile", "momfile", "jfile", "inpsd.minimal")


def rspt_get_input_data(args):
    """
    Read the input data for RSPt calculations to get structure info.

    Author:
        Anders Bergman
    """
    types, pos = extract_position_data(args.data)
    return types, pos


def rspt_relativistic_workflow(args):
    """
    Main function to extract and process relativistic exchange data
    from RSPT output files.

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

    print(50 * "-")
    print("Relativistic exchange workflow")
    # Define input directions and atoms
    input_directions = ["100", "010", "001"]
    atomlist = args.atoms
    input_files = {
        input_directions[0]: args.rx,
        input_directions[1]: args.ry,
        input_directions[2]: args.rz,
    }

    exchange_data: dict[str, list[RsptExchange]] = {}

    # Extract information from each input file
    for direction in input_directions:
        print("Extracting exchange data from:", direction)
        exchange_data[direction] = []
        for file_name in input_files[direction]:
            # print(f"Atom: {atoms}")
            extracted_data = RsptExchange(file_name, args.maptype)
            exchange_data[direction].append(extracted_data.outmap)

    concatenated_data: dict[str, np.ndarray] = {}
    for key, value in exchange_data.items():
        concatenated_data[key] = np.concatenate(value)

    input_atoms = np.unique(concatenated_data["001"][:, 0]).astype(np.int32).tolist()

    # Find the indices of input_atoms that match the entries in atomlist
    print(f"Input atoms: \n  {input_atoms}")
    if atomlist is not None:
        output_atoms = [i for i, atom in enumerate(input_atoms) if atom in atomlist]
    else:
        # output_atoms = [ i + 1 for i in input_atoms]
        output_atoms = list(range(len(input_atoms)))
    print(f"Output atoms: \n  {output_atoms}")

    # Filter data according to available atoms
    if atomlist:
        print(f"Atomlist: {atomlist}")
        filter_list = [i - 1 for i in atomlist]
    else:
        filter_list = [i - 1 for i in input_atoms]

    print(f"Filtered atoms: \n  {filter_list}")

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

    # Exctract data from an scf-run
    print(50 * "-")
    print("Extracting SCF data")
    if args.scf and args.data:
        scf_data = RsptData(args)
        # Print the lattice, basis, and moments information
        scf_data.print_lattice("lattice.dat")
        scf_data.print_positions("posfile")
        scf_data.print_moments("momfile")
        scf_data.print_template("posfile", "momfile", "jfile", "inpsd.minimal")
    else:
        print(
            "No SCF and/or `data` file provided, full `UppASD` template not generated."
        )

    # Filter out the Jij values so that only selected atoms are included
    # eg. only the magnetic atoms
    # j_truncated = downscale_exchange(summed_data, [1, 2])
    if atomlist is None:
        # j_truncated = downscale_exchange(concatenated_data, input_atoms)
        print("Thresholded atoms:", scf_data.filter_list)
        j_filtered = downscale_exchange(summed_data, scf_data.filter_list)
        # j_truncated = np.array(summed_data, dtype=np.float64)
    else:
        j_filtered = downscale_exchange(summed_data, atomlist)

    # Create the various projections you can make from the exchange interactions
    j_dict = extract_projections(j_filtered)
    # Print all projections to file
    print_projections(j_dict, args.maptype, is_relativistic=True)

    print("Extraction and storage completed successfully!")
    print(50 * "-")


def rspt_scalar_workflow(args):
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
    # rspt_exchange = RsptExchange(args)
    # rspt_exchange.save_exchange()
    print(50 * "-")
    print("Scalar exchange workflow")

    # Process the exchange interactions
    input_files = args.exchange
    atomlist = args.atoms
    exchange_list: list[RsptExchange] = []
    exchange_data = []

    for file_name in input_files:
        extracted_data = RsptExchange(file_name, args.maptype)
        exchange_list.append(extracted_data)
        exchange_data.append(extracted_data.outmap)

    concatenated_data = np.concatenate(exchange_data)
    input_atoms = np.unique(concatenated_data[:, 0]).astype(np.int32).tolist()

    # Find the indices of input_atoms that match the entries in atomlist
    print("Input atoms:", input_atoms)
    if atomlist is not None:
        output_atoms = [i for i, atom in enumerate(input_atoms) if atom in atomlist]
    else:
        # output_atoms = [ i + 1 for i in input_atoms]
        output_atoms = list(range(len(input_atoms)))
    print("Output atoms:", output_atoms)

    # Filter data according to available atoms
    if atomlist:
        print(f"Atomlist: {atomlist}")
        filter_list = [i - 1 for i in atomlist]
    else:
        filter_list = [i - 1 for i in input_atoms]

    print(f"Filtered atoms: \n{filter_list}")

    # Exctract data from an scf-run
    print("-" * 50)
    if args.scf and args.data:
        scf_data = RsptData(args)
        # Print the lattice, basis, and moments information
        scf_data.print_lattice("lattice.dat")
        scf_data.print_positions("posfile")
        scf_data.print_moments("momfile")
        scf_data.print_template("posfile", "momfile", "jfile", "inpsd.minimal")
        # lattice = scf_data.lattice
        # basis = scf_data.basis
        # moments = scf_data.moments
        # atoms = scf_data.atoms
    else:
        # lattice = exchange_list[0].lattice
        # basis = [obj.r_i for obj in exchange_list]
        # moments = [1.0 for obj in exchange_list]
        print(
            "No SCF and/or `data` file provided, full `UppASD` template not generated."
        )

    # Filter out the Jij values so that only selected atoms are included
    # eg. only the magnetic atoms
    if atomlist is None:
        # j_truncated = downscale_exchange(concatenated_data, input_atoms)
        print("Thresholded atoms:", scf_data.filter_list)
        j_filtered = downscale_exchange(concatenated_data, scf_data.filter_list)
    else:
        print("Atomlist:", atomlist)
        j_filtered = downscale_exchange(concatenated_data, atomlist)

    # Exctact all projections from the exchange interactions
    j_extracted = extract_projections(j_filtered, is_relativistic=False)

    # Print all projections to file
    print_projections(
        j_extracted, maptype=args.maptype, is_relativistic=False, cutoff=args.cutoff
    )

    print("Extraction and storage completed successfully!")
    print(50 * "-")


def main():
    """
    Main function to parse command-line arguments and run the appropriate
    workflow based on the provided options.

    This function sets up a command-line interface (CLI) with three options:
    - `-s` or `--scf`: Run the SCF workflow with the given file path.
    - `-e` or `--exchange`: Run the exchange workflow with the given file path.
    - `-d` or `--dir`: Run the default workflow with the given directories.

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
    parser.add_argument("-d", "--data", type=str, help="FILE_NAME for RSPt input data")
    parser.add_argument(
        "-s", "--scf", type=str, help="FILE_NAME for RSPt SCF output data"
    )
    parser.add_argument(
        "-e",
        "--exchange",
        type=str,
        nargs="+",
        metavar="FILE_NAME",
        help="Run exchange workflow for the given FILE_NAMEs",
    )
    parser.add_argument(
        "--rx",
        type=str,
        nargs="+",
        metavar="FILE_NAME",
        help="Input files for relativistic exchange workflow (x-direction)",
    )
    parser.add_argument(
        "--ry",
        type=str,
        nargs="+",
        metavar="FILE_NAME",
        help="Input files for relativistic exchange workflow (y-direction)",
    )
    parser.add_argument(
        "--rz",
        type=str,
        nargs="+",
        metavar="FILE_NAME",
        help="Input files for relativistic exchange workflow (z-direction)",
    )
    parser.add_argument(
        "-m",
        "--maptype",
        type=str,
        default="C",
        help="Output maptype: (C)artesian, (D)irect or maptype (3)",
    )
    parser.add_argument(
        "-a",
        "--atoms",
        nargs="+",
        type=int,
        help="List of atoms to extract exchange interactions for",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.0,
        help="Threshold for moment magnitudes. Default is 0.0",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        help="Radius for excange interaction cutoff",
    )

    # Parse the arguments
    args = parser.parse_args()

    print(50 * "-")
    print(10 * " " + "RSPt Exchange Extraction Tool")
    print(50 * "-")
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    # Decide which function to run based on the CLI arguments
    if args.exchange:
        rspt_scalar_workflow(args)
    elif args.rx and args.ry and args.rz:
        rspt_relativistic_workflow(args)
    elif args.scf:
        rspt_get_scf_data(args)


if __name__ == "__main__":
    main()
