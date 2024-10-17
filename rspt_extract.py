"""
This module contains core routines for extracting data from RPSt output.

The functions in this module are designed to extract various data from files generated
by RPSt calculations.  These calculations involve extracting and transforming matrices,
extracting lattice constants and matrices, extracting basis vectors, extracting distance
vectors, extracting exchange matrices, and converting atomic coordinates to maptype three.

The main functions in this module are:
- flatten_and_concatenate: Flattens and concatenates elements of a list.
- transform_matrix: Transforms a given matrix according to a specific formula.
- extract_bravais_lattice_matrix: Extracts the bravais lattice matrix from a file.
- extract_basis_vectors: Extracts basis vectors from a file.
- extract_distance_vectors: Extracts vectors from a file.
- extract_exchange_matrices: Extracts matrices from a file containing calculated quantities.
- convert_to_maptype_three: Converts atomic coordinates to maptype three.

Each function is documented with its parameters, return types, and a
brief description of its functionality.
"""

import numpy as np


def flatten_and_concatenate(input_list):
    """
    Flattens and concatenates the elements of the input list.

    Args:
        input_list (list): The input list containing elements to be flattened and concatenated.

    Returns:
        numpy.ndarray: The flattened and concatenated array.
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
    Transforms the given matrix according to a specific formula.

    Args:
        matrix (numpy.ndarray): The input matrix.

    Returns:
        numpy.ndarray: The transformed matrix.
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


def extract_bravais_lattice_matrix(file_path):
    """
    Extracts the bravais lattice matrix from the given file.

    Args:
        file_path (str): The path to the file.

    Returns:
        tuple: A tuple containing the lattice constant and the lattice matrix.
    """
    print("Extracting bravais lattice matrix from:", file_path)
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    matrix = []
    capture = False
    alat = 0.0

    for line in lines:
        if "Bravais lattice basis" in line:
            capture = True
            alat = np.float32(line.split()[4])
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


def extract_moments(file_path):
    """
    Extracts magnetic moments from a file.

    Parameters:
    file_path (str): The path to the file containing the basis vectors.

    Returns:
    numpy.ndarray: An array of magnetic moments

    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    moments = []
    capture = False
    captured = False

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

    return np.array(moments, dtype=np.float32)


def extract_basis_vectors(file_path):
    """
    Extracts basis vectors from a file.

    Parameters:
    file_path (str): The path to the file containing the basis vectors.

    Returns:
    numpy.ndarray: An array of basis vectors.

    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    vectors = []
    capture = False

    for line in lines:
        if "G*t/2pi" in line:
            capture = True
            continue  # Skip the line containing 'Bravais lattice basis'
        if capture:
            # Check if the line has three floating point numbers
            values = line.split()
            vectors.append([float(val) for val in values[0:3]])
            capture = False

    return np.array(vectors, dtype=np.float32)


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

    """
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    r_ij = []
    d_ij = []
    i_atom = 0
    j_atom = []
    r_ij = []
    capture_glob = False

    for line in lines:
        if "Central" in line:
            capture_glob = True
            i_atom = np.int32(line.split()[1])
            r_i = np.array(line.split()[5:8], dtype=np.float32)
            continue

        if "END" in line:
            capture_glob = False

        if capture_glob:
            values = line.split()
            r_ij.append(np.array(values[5:8], dtype=np.float32))
            d_ij.append(np.float32(values[10]))
            j_atom.append(np.int32(values[1]))

    return i_atom, r_i, j_atom, d_ij, r_ij


def extract_exchange_matrices(file_path):
    """
    Extracts matrices from a file containing calculated quantities in a global coordinate system.

    Args:
        file_path (str): The path to the file containing the calculated quantities.

    Returns:
        list: A list of numpy arrays representing the extracted matrices.

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


def convert_to_maptype_three(i_atom, j_atoms, basis, lattice, alat, r_ij):
    """
    Converts the given atomic coordinates to maptype three.

    Parameters:
    i_atom (int): Index of the central atom.
    j_atoms (list): List of indices of the neighboring atoms.
    basis (list): List of atomic coordinates.
    lattice (ndarray): Lattice vectors.
    alat (float): Lattice constant.
    r_ij (ndarray): Relative positions of the neighboring atoms.

    Returns:
    ndarray: Array of converted atomic coordinates.

    """
    invlatt = np.linalg.inv(lattice)
    basis[:, 2] = -basis[:, 2]
    s_ij = []
    for idx, vector in enumerate(r_ij):
        r_i = np.dot(lattice, basis[i_atom - 1])
        r_j = np.dot(lattice, basis[j_atoms[idx] - 1])
        r_shift = vector / alat
        v_ij = r_shift - r_j  # - r_i
        s_ij.append(np.dot(invlatt, v_ij))

    return np.array(s_ij, dtype=np.float32).round(4)
