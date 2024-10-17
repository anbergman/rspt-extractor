"""
This module provides functionality to extract, process, and save exchange data 
from a specified file. It includes the RsptExchange class and several utility 
functions for handling exchange matrices and projections.
Classes:
    RsptExchange: A class to handle extraction, processing, and saving of 
    exchange data.
Functions:
    downscale_exchange(exchange: list, mask_list: list) -> np.ndarray:
    extract_projections(exchange: np.ndarray) -> dict:
        Extract various projections from the exchange matrix.
    print_projections(j_dict: dict) -> None:

Author:
    Anders Bergman
"""
import numpy as np
from .rspt_extract import (
    extract_bravais_lattice_matrix,
    extract_exchange_matrices,
    extract_distance_vectors,
    extract_basis_vectors,
    convert_to_maptype_three,
    transform_matrix,
    flatten_and_concatenate,
)
# import .rspt_extract as rs


class RsptExchange:
    """
        This class provides methods to extract exchange data from a specified file,
        process the data, and save the processed data to an output file.

        Attributes:
            alat (float): Lattice constant.
            lattice (np.ndarray): Bravais lattice matrix.
            matrices (np.ndarray): Exchange matrices.
            i_atom (np.ndarray): Atom indices.
            r_i (np.ndarray): Position vectors of the i-th atom.
            j_atoms (np.ndarray): Atom indices of neighboring atoms.
            d_ij (np.ndarray): Distance vectors between atoms.
            r_ij (np.ndarray): Distance vectors between atoms in reduced coordinates.
            basis (np.ndarray): Basis vectors.
            s_ij (np.ndarray): Maptype three representation of distance vectors.
            outmap (np.ndarray): Processed data.

        Methods:
            __init__(file_path: str):
                Initializes the RsptExchange object with the given file path.
            
            extract_data():
                Extracts exchange data from the input file.
            
            process_data():
                Processes the extracted data.
            
            save_output(output_path: str):
                Saves the processed data to a specified output file.

        Example:
            file_path = "spin-001/out-3"
            output_path = "outmap.txt"

            rspt_exchange = RsptExchange(file_path)
            rspt_exchange.extract_data()
            rspt_exchange.process_data()
            rspt_exchange.save_output(output_path)

        Author:
            Anders Bergman
    """

    def __init__(self, file_path):
        """
            This constructor sets up the RsptExchange object by initializing its 
            attributes and calling methods to extract and process data from the 
            provided input file.


            Attributes:
                alat (float or None): Lattice parameter.
                lattice (list or None): Lattice vectors.
                matrices (list or None): Exchange matrices.
                i_atom (int or None): Index of the atom.
                r_i (list or None): Position vector of the atom.
                j_atoms (list or None): Indices of neighboring atoms.
                d_ij (list or None): Distance vectors between atoms.
                r_ij (list or None): Relative position vectors.
                basis (list or None): Basis vectors.
                s_ij (list or None): Spin interaction parameters.
                outmap (dict or None): Output mapping of processed data.

            Author: Anders Bergman
        """
        self.file_path = file_path
        self.alat = None
        self.lattice = None
        self.matrices = None
        self.i_atom = None
        self.r_i = None
        self.j_atoms = None
        self.d_ij = None
        self.r_ij = None
        self.basis = None
        self.s_ij = None
        self.outmap = None
        self.extract_data()
        self.process_data()

    def extract_data(self):
        """
        Extracts various data from the specified file and assigns it to instance 
        variables.

        This method performs the following operations:
        1. Extracts the Bravais lattice matrix and lattice constant.
        2. Extracts the exchange matrices.
        3. Extracts distance vectors including atom indices, positions, and 
           distances.
        4. Extracts basis vectors.
        5. Converts the extracted data to a specific map type.

        The extracted data is stored in the following instance variables:
        - alat: Lattice constant.
        - lattice: Bravais lattice matrix.
        - matrices: Exchange matrices.
        - i_atom: Indices of atoms.
        - r_i: Positions of atoms.
        - j_atoms: Indices of neighboring atoms.
        - d_ij: Distance vectors between atoms.
        - r_ij: Distance magnitudes between atoms.
        - basis: Basis vectors.
        - s_ij: Converted map type data.

        Args:
            None

        Returns:
            None

        Author: Anders Bergman
        """
        self.alat, self.lattice = extract_bravais_lattice_matrix(
            self.file_path)
        self.matrices = extract_exchange_matrices(self.file_path)
        self.i_atom, self.r_i, self.j_atoms, self.d_ij, self.r_ij = (
            extract_distance_vectors(self.file_path))
        self.basis = extract_basis_vectors(self.file_path)
        self.s_ij = convert_to_maptype_three(self.i_atom, self.j_atoms,
                                             self.basis, self.lattice,
                                             self.alat, self.r_ij)

    def process_data(self):
        """
        Processes the data to generate a flattened and concatenated output map.

        This method iterates over the `r_ij` attribute, constructs a list of 
        transformed and normalized data for each index, and appends the result 
        to the `outmap` attribute. The final `outmap` is stored as a NumPy array.

        Args:
            None

        Returns:
            None

        Attributes:
            outmap (np.ndarray): A NumPy array containing the processed data.

        Raises:
            None

        Author:
            Anders Bergman
        """
        outmap = []
        for idx, _ in enumerate(self.r_ij):
            i_list = [
                self.i_atom,
                self.j_atoms[idx],
                self.s_ij[idx],
                transform_matrix(self.matrices[idx]).flatten(),
                self.d_ij[idx] / self.alat,
            ]
            outmap.append(flatten_and_concatenate(i_list))
        self.outmap = np.array(outmap)

    def save_output(self, output_path):
        """
        Save the output map to a specified file.

        This method saves the `self.outmap` array to the given `output_path` using
        a specific format for the data. The format includes both integer and float
        values with defined precision.

        Args:
            output_path (str): The path to the file where the output will be saved.

        Returns:
            None

        Example:
            >>> obj.save_output('/path/to/output/file.txt')

        Author:
            Anders Bergman
        """
        fmt1 = "%4d %4d   % 4.1f % 4.1f % 4.1f   "
        fmt2 = "% 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f     %8.4f"
        fmt = fmt1 + fmt2
        np.savetxt(output_path, self.outmap, fmt=fmt)


def downscale_exchange(exchange, mask_list):
    """
    Downscale the exchange matrix based on the provided mask list.

    This function filters and modifies the given exchange matrix such that only
    the rows where both elements are present in the mask list are retained.
    The indices of these elements are then updated to reflect their positions
    in the mask list.

    Args:
        exchange (list of list of int): The exchange matrix to be downscaled.
            Each row is expected to be a list of integers where the first two
            elements represent indices.
        mask_list (list of int): The list of indices to be used for masking
            and downscaling the exchange matrix.

    Returns:
        np.ndarray: A numpy array of the downscaled exchange matrix with dtype
        np.float32.

    Author:
        Anders Bergman
    """
    masked_data = []
    for row in exchange:
        if row[0] in mask_list and row[1] in mask_list:
            new_row = row.copy()
            new_row[0] = mask_list.index(row[0]) + 1
            new_row[1] = mask_list.index(row[1]) + 1
            masked_data.append(new_row)

    return np.array(masked_data, dtype=np.float32)


def extract_projections(exchange):
    """
    Extracts various projections from the exchange matrix.

    This function processes the exchange matrix to extract scalar, diagonal, 
    tensor, Dzyaloshinskii-Moriya interaction (DMI), and symmetric components 
    for each atom. It also separates the left and right parts of the exchange 
    matrix.

    Args:
        exchange (np.ndarray): A 2D numpy array representing the exchange 
            matrix. The shape of the array is (natom, nfeatures).

    Returns:
        dict: A dictionary containing the following keys:
            - "scalar" (np.ndarray): Scalar projection for each atom.
            - "diagonal" (np.ndarray): Diagonal projection for each atom.
            - "tensor" (np.ndarray): Tensor projection for each atom.
            - "dmi" (np.ndarray): DMI projection for each atom.
            - "symmetric" (np.ndarray): Symmetric projection for each atom.
            - "left" (np.ndarray): Left part of the exchange matrix.
            - "right" (np.ndarray): Right part of the exchange matrix.

    Author: Anders Bergman
    """
    natom = exchange.shape[0]
    j_dict = {}
    j_dict["scalar"] = np.zeros((natom, 1))
    j_dict["diagonal"] = np.zeros((natom, 9))
    j_dict["tensor"] = np.zeros((natom, 9))
    j_dict["dmi"] = np.zeros((natom, 3))
    j_dict["symmetric"] = np.zeros((natom, 3))
    for irow, row in enumerate(exchange):
        j_mat = row[5:14].reshape(3, 3)
        j_mat_sym = 0.5 * (j_mat - j_mat.T)
        j_mat_asym = 0.5 * (j_mat - j_mat.T)
        dmi = np.array([j_mat_asym[1, 2], j_mat_asym[2, 0], j_mat_asym[0, 1]])
        jsy = np.array([j_mat_sym[1, 2], j_mat_sym[2, 0], j_mat_sym[0, 1]])
        j_dict["scalar"][irow] = np.trace(j_mat) / 3.0
        j_dict["diagonal"][irow] = np.diag(np.diag(j_mat)).flatten()
        j_dict["dmi"][irow] = dmi
        j_dict["symmetric"][irow] = jsy
        j_dict["tensor"][irow] = j_mat.flatten()

    j_dict["left"] = exchange[:, 0:5]
    j_dict["right"] = exchange[:, 14:]

    return j_dict


def print_projections(j_dict):
    """
    Save projections from a dictionary to files with specific formatting.

    This function takes a dictionary containing projection data and saves the 
    projections to files. Each type of projection is saved in a separate file 
    with a specific format.

    Args:
        j_dict (dict): A dictionary containing the projection data. The keys 
        should include "left", "right", and one or more of the following: 
        "scalar", "diagonal", "dmi", "symmetric", "tensor". Each key should 
        map to a numpy array.

    Returns:
        None

    Files Created:
        - j_scalar.dat
        - j_diagonal.dat
        - j_dmi.dat
        - j_symmetric.dat
        - j_tensor.dat

    Example:
        j_dict = {
            "left": np.array([...]),
            "right": np.array([...]),
            "scalar": np.array([...]),
            ...
        }
        print_projections(j_dict)

    Author:
        Anders Bergman
    """
    fmt_l = "%4d %4d   % 4.1f % 4.1f % 4.1f    "
    fmt_1 = "% 10.6f"
    fmt_3 = "% 10.6f % 10.6f % 10.6f"
    fmt_9 = "% 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f"
    fmt_r = "     %8.4f"

    keys = ["scalar", "diagonal", "dmi", "symmetric", "tensor"]
    fmts = [fmt_1, fmt_9, fmt_3, fmt_3, fmt_9]

    for idx, key in enumerate(keys):
        outmat = np.hstack((j_dict["left"], j_dict[key], j_dict["right"]))
        fname = f"j_{key}.dat"
        np.savetxt(fname, outmat, fmt=fmt_l + fmts[idx] + fmt_r)
    print("Projections saved to files.")
    return


if __name__ == "__main__":
    # Example usage
    FILE_PATH = "spin-001/out-3"
    OUTPUT_PATH = "outmap.txt"

    rspt_exchange = RsptExchange(FILE_PATH)
    rspt_exchange.save_output(OUTPUT_PATH)
