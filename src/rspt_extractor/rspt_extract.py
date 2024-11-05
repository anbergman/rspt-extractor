"""
Module for extracting and transforming data from RSPT output files.

This module provides functions to extract various types of data from RSPT
output files, including Bravais lattice matrices, magnetic moments, basis
vectors, distance vectors, and exchange matrices. It also includes functions
to transform matrices and convert atomic positions to a specific map type.

Functions:
    flatten_and_concatenate(input_list: list) -> np.ndarray:

    transform_matrix(matrix: np.ndarray) -> np.ndarray:

    extract_lattice_scf(file_path: str) -> tuple:

    extract_moments_scf(file_path: str) -> np.ndarray:

    extract_basis_scf(file_path: str) -> np.ndarray:

    extract_distance_vectors(file_path: str) -> tuple:

    extract_exchange_matrices(file_path: str) -> list:

    convert_to_maptype_three(i_atom, j_atoms: list,
                            basis: np.ndarray, lattice: np.ndarray,
                            alat: float, r_ij: np.ndarray) -> np.ndarray:
"""

import re
import numpy as np


def flatten_and_concatenate(input_list):
    """
    Flattens and concatenates elements from the input list.

    This function takes a list of elements, where each element can be either
    a numpy ndarray or any other type. If an element is an ndarray, it is
    flattened and its elements are added to the result list. If an element is
    not an ndarray, it is directly added to the result list. The final result
    is returned as a numpy array.

    Args:
        input_list (list): A list containing elements that are either numpy
            ndarrays or other types.

    Returns:
        np.ndarray: A numpy array containing the flattened and concatenated
            elements from the input list.

    Example:
        >>> import numpy as np
        >>> input_list = [np.array([1, 2]), 3, np.array([4, 5])]
        >>> flatten_and_concatenate(input_list)
        array([1, 2, 3, 4, 5])

    Author:
        Anders Bergman
    """
    flattened_elements = []
    for item in input_list:
        if isinstance(item, np.ndarray):
            flattened_elements.extend(item.flatten())
        else:
            flattened_elements.append(item)
    return np.array(flattened_elements)


def transform_matrix(matrix):
    """
    Transforms a 3x3 matrix according to a specific formula.

    This function takes a 3x3 matrix and applies a transformation to it,
    resulting in a new 3x3 matrix. The transformation is defined by the
    following formula:

    transformed_matrix = [
        [jx, cz + dz, cy - dy],
        [cz - dz, jy, cx + dx],
        [cy + dy, cx - dx, jz]
    ]

    Args:
        matrix (np.ndarray): A 3x3 numpy array representing the original matrix.

    Returns:
        np.ndarray: A 3x3 numpy array representing the transformed matrix.

    Example:
        >>> import numpy as np
        >>> matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> transform_matrix(matrix)
        array([[ 1, 15, -1],
               [-1,  5,  5],
               [11, -1,  9]])

    Author:
        Anders Bergman
    """
    # Extract the elements from the original matrix
    jx, dx, cx = matrix[0]
    jy, dy, cy = matrix[1]
    jz, dz, cz = matrix[2]

    # Create the transformed matrix according to the given formula
    transformed_matrix = np.array(
        [[jx, cz + dz, cy - dy], [cz - dz, jy, cx + dx], [cy + dy, cx - dx, jz]]
    )

    return transformed_matrix


def check_relativistic(file_name):
    """
    Check if the calculation is relativistic.

    This function reads a file and checks if the calculation is relativistic
    by searching for the keyword "fulrel" in the file.

    Args:
        file_name (str): The path to the file containing the calculation data.

    Returns:
        bool: True if the calculation is relativistic, False otherwise.

    Raises:
        FileNotFoundError: If the file specified by `file_name` does not exist.

    Example:
        >>> is_relativistic = check_relativistic("path/to/file.txt")
        >>> print(is_relativistic)
    """
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if "fulrel" in line:
            return lines[i + 1].strip().split()[0] == "T"

    return False


def extract_lattice_scf(file_name):
    """
    Extracts the Bravais lattice matrix from a given file.

    This function reads a file specified by `file_path` and extracts the Bravais
    lattice matrix along with the lattice constant (alat). The file is expected
    to contain a section labeled "Bravais lattice basis" followed by three lines
    of floating point numbers representing the matrix.

    Args:
        file_path (str): The path to the file containing the Bravais lattice
                         matrix data.

    Returns:
        tuple: A tuple containing:
            - alat (float): The lattice constant.
            - matrix (numpy.ndarray): A 3x3 numpy array representing the Bravais
                                      lattice matrix.

    Raises:
        FileNotFoundError: If the file specified by `file_path` does not exist.
        ValueError: If the file does not contain the expected Bravais lattice
                    matrix data.

    Example:
        >>> alat, matrix = extract_bravais_lattice_matrix("path/to/file.txt")
        >>> print(alat)
        >>> print(matrix)

    Author:
        Anders Bergman
    """
    # print("Extracting bravais lattice matrix from:", file_path)
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    matrix = []
    capture = False
    alat = 0.0

    for line in lines:
        if "Bravais lattice basis" in line:
            capture = True
            alat = np.float64(line.split()[4])
            continue  # Skip the line containing 'Bravais lattice basis'
        if capture:
            # Check if the line has three floating point numbers
            values = line.split()
            if len(values) == 3:
                matrix.append([float(val) for val in values])

            # Stop capturing after reading three lines
            if len(matrix) == 3:
                break

    return alat, np.array(matrix)


def extract_moments_scf(file_name, is_relativistic):
    """
    Extracts magnetic moments from a specified file.

    This function reads a file containing magnetic moment data and extracts the
    spin, orbital, and total moments. The data is expected to be in a specific
    format where the line containing
    "type           spin     orbital       total"
    marks the beginning of the relevant data section for relativistic data, and
    "type        spin"
    for non-relativistic data.

    Args:
        file_path (str): The path to the file containing the magnetic moment data.

    Returns:
        np.ndarray: A NumPy array of shape (n, 3) where n is the number of data
        points. Each row contains the spin, orbital, and total moments as
        floating-point numbers.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the data format in the file is incorrect.

    Example:
        >>> moments = extract_moments("path/to/datafile.txt")
        >>> print(moments)
        [[0.1, 0.2, 0.3],
         [0.4, 0.5, 0.6],
         ...]

    Author:
        Anders Bergman
    """
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    moments = []
    capture = False
    captured = False

    if is_relativistic:
        for line in lines:
            if "type           spin     orbital       total" in line:
                capture = True
                continue

            if capture and not captured:
                # Check if the line has three floating point numbers
                values = line.split()
                if values[0] == "total":
                    captured = True
                    capture = False
                    continue
                moments.append([float(val) for val in values[1:4]])
    else:
        for line in lines:
            if "type        spin" in line:
                capture = True
                continue

            if capture and not captured:
                # Check if the line has three floating point numbers
                values = line.split()
                if len(values) == 0:
                    captured = True
                    capture = False
                    continue
                moments.append([float(values[1]), 0.0, float(values[1])])

    return np.array(moments, dtype=np.float64)


def extract_basis_scf(file_name):
    """
    Extracts basis vectors from a file.

    Parameters:
    file_path (str): The path to the file containing the basis vectors.

    Returns:
    numpy.ndarray: An array of basis vectors.

    """
    with open(file_name, "r", encoding="utf-8") as file:
        lines = file.readlines()

    vectors = []
    capture = False

    # Define a regular expression to match floats (including scientific notation)
    float_pattern = r"([+-]?\d+\.\d+E[+-]?\d+)"

    for line in lines:
        if "G*t/2pi" in line:
            capture = True
            continue  # Skip the line containing 'Bravais lattice basis'
        if capture:
            # Check if the line has three floating point numbers
            # values = line.split()

            # Use re.findall to extract all the float numbers from the string
            parsed_data = re.findall(float_pattern, line)

            # Convert strings to float (optional)
            parsed_data = [float(num) for num in parsed_data]

            vectors.append([float(val) for val in parsed_data])
            # vectors.append([float(val) for val in values[0:3]])
            capture = False

    return np.array(vectors, dtype=np.float64)


def extract_distance_vectors(file_path):
    """
    Extracts vectors from a file.

    Args:
        file_path (str): The path to the file.

    Returns:
        tuple: A tuple containing the following elements:
            - i_atom (numpy.int32): The central atom index.
            - r_i (numpy.ndarray): The coordinates of the central atom.
            - j_atom (list): A list of atom indices.
            - d_ij (list): A list of distances.
            - r_ij (list): A list of coordinates.

    Note:
        The distance vectors r_ij, and the coordinates r_i are in the global
        cartesian coordinate system, scaled by the lattice constant.

    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    r_ij = []
    d_ij = []
    i_atom = 0
    r_i = 0.0
    j_atom = []
    r_ij = []
    capture_glob = False

    for line in lines:
        if "Central" in line:
            capture_glob = True
            i_atom = np.int32(line.split()[1])
            r_i = np.array(line.split()[5:8], dtype=np.float64)
            continue

        if "END" in line:
            capture_glob = False

        if capture_glob:
            values = line.split()
            r_ij.append(np.array(values[5:8], dtype=np.float64) - r_i)
            d_ij.append(np.float64(values[10]))
            j_atom.append(np.int32(values[1]))

    return i_atom, r_i, j_atom, d_ij, r_ij


def extract_position_data(filename):
    """
    Reads RSPt data from a given file and extracts atom types and positions.

    This function parses a file to extract the number of atoms and their
    positions. It looks for specific keywords ('natom' and 'tau1') to identify
    the relevant sections of the file. The positions are extracted as a list
    of coordinates, and the types are assigned based on the order of appearance.

    Args:
        filename (str): The path to the file containing RSPT data.

    Returns:
        Tuple[List[int], List[List[float]]]: A tuple containing two lists:
            - types: A list of integers representing the type of each atom.
            - positions: A list of lists, where each sublist contains three
              floats representing the x, y, and z coordinates of an atom.

    Raises:
        ValueError: If the file format is incorrect or if 'natom' is not found
                    before 'tau1'.

    Author:
        Anders Bergman
    """
    types = []
    positions = []
    ntype = 0
    with open(filename, "r", encoding="utf-8") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            # Look for the line containing 'natom'
            if "natom" in line:
                # The value is in the line below, first column
                natom = int(lines[i + 1].split()[0])
                ntype += 1
            elif "tau1" in line:
                for iatom in range(natom):
                    # The positions are in the lines below,
                    # starting from the line after 'tau1'
                    positions.append(
                        [float(val) for val in lines[i + 1 + iatom].split()[0:3]]
                    )
                    types.append(ntype)

    return types, positions


def extract_exchange_matrices(file_path):
    """
    Extracts exchange matrices from a specified file.

    This function reads a file containing calculated quantities in the global
    coordinate system and extracts exchange matrices. The matrices are
    identified by specific markers in the file and are returned as a list of
    NumPy arrays.

    Args:
        file_path (str): The path to the file containing the exchange matrices.

    Returns:
        List[np.ndarray]: A list of 3x3 NumPy arrays representing the exchange
        matrices.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file contains invalid data that cannot be parsed
        into matrices.

    Example:
        >>> matrices = extract_exchange_matrices("path/to/file.txt")
        >>> for matrix in matrices:
        >>>     print(matrix)

    Author:
        Anders Bergman
    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    matrices = []
    capture_glob = False
    capture_mats = False
    current_matrix = []

    for line in lines:
        if "CALCULATED QUANTITIES IN GLOBAL COORDINATE SYSTEM" in line:
            capture_glob = True
        if "END" in line:
            capture_glob = False
            capture_mats = False
        if "Neighbor   dist." in line:
            capture_mats = True
            continue

        if capture_glob and capture_mats:
            if line.strip().startswith(("x", "y", "z")):
                values = line.split()
                if len(values) == 4:
                    _, j_ij, d_ij, c_ij = values
                    current_matrix.append([float(j_ij), float(d_ij), float(c_ij)])
                    if len(current_matrix) == 3:
                        matrices.append(np.array(current_matrix))
                        current_matrix = []

    return matrices


def extract_exchange_scalars(file_path):
    """
    Extracts exchange interactions from a specified file.

    This function reads a file containing calculated quantities in the global
    coordinate system and extracts exchange constants. The constants are
    identified by specific markers in the file and are returned as a list of
    NumPy arrays.

    Args:
        file_path (str): The path to the file containing the exchange matrices.

    Returns:
        List[np.ndarray]: A list of 3x3 NumPy arrays representing the exchange
        matrices.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file contains invalid data that cannot be parsed
        into matrices.

    Example:
        >>> matrices = extract_exchange_matrices("path/to/file.txt")
        >>> for matrix in matrices:
        >>>     print(matrix)

    Author:
        Anders Bergman
    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    couplings = []
    capture = False

    for line in lines:
        if "Neighbor     dist." in line:
            capture = True
            continue
        if "END" in line:
            capture = False

        if capture:
            couplings.append(float(line.split()[-1]))

    return couplings


def convert_to_maptype_three(i_atom, j_atoms, basis, lattice, alat, r_ij):
    """
    Converts atomic positions to a specific map type (type three).

    This function takes in atomic positions and lattice parameters to compute
    the mapped positions in a new coordinate system.

    Args:
        j_atoms (list[int]): List of atom indices.
        basis (np.ndarray): Basis vectors of the atomic positions.
        lattice (np.ndarray): Lattice vectors of the crystal.
        alat (float): Lattice constant.
        r_ij (np.ndarray): Relative positions between atoms.

    Returns:
        np.ndarray: Transformed atomic positions in the new coordinate system,
        rounded to 4 decimal places.

    Author:
        Anders Bergman
    """
    invlatt = np.linalg.inv(lattice)
    # basis[:, 2] = -basis[:, 2]
    s_ij = []
    for idx, vector in enumerate(r_ij):
        r_i = np.dot(lattice, basis[i_atom - 1])
        r_j = np.dot(lattice, basis[j_atoms[idx] - 1])
        r_shift = vector / alat
        v_ij = r_shift - (r_j - r_i)
        s_ij.append(np.dot(invlatt, v_ij))

    return np.array(s_ij, dtype=np.float64).round(4)


def convert_to_direct(lattice, alat, r_ij):
    """
    Converts cartesian atomic positions direct coordinates

    This function takes in atomic positions and lattice parameters to compute
    the mapped positions in a new coordinate system.

    Args:
        lattice (np.ndarray): Lattice vectors of the crystal.
        alat (float): Lattice constant.
        r_ij (np.ndarray): Relative positions between atoms.

    Returns:
        np.ndarray: Transformed atomic positions in the new coordinate system,
        rounded to 4 decimal places.

    Author:
        Anders Bergman
    """
    invlatt = np.linalg.inv(lattice)
    # basis[:, 2] = -basis[:, 2]
    s_ij = []
    for vector in r_ij:
        r_shift = vector / alat
        s_ij.append(np.dot(invlatt.T, r_shift))

    return np.array(s_ij, dtype=np.float64).round(4)
