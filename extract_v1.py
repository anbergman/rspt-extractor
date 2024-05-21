import numpy as np

def flatten_and_concatenate(input_list):
    flattened_elements = []
    for item in input_list:
        if isinstance(item, np.ndarray):
            flattened_elements.extend(item.flatten())
        else:
            flattened_elements.append(item)
    return np.array(flattened_elements)

def modshift(arr):
    """
    Perform modulo operation on a (3,) ndarray to ensure all entries are between 0.0 and 1.0

    Parameters:
        arr (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The modified array.
    """
    return np.mod(arr, 1)

def modshift_r(arr):
    """
    Perform modulo operation on a (3,) ndarray to ensure all entries are between -0.5 and +0.5.

    Parameters:
        arr (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The modified array.
    """
    return np.mod(arr + 0.5, 1) - 0.5

def transform_matrix(matrix):
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
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    matrix = []
    capture = False

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

def extract_basis_vectors(file_path):
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

    return np.array(vectors,dtype = np.float32)


def extract_vectors(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    r_ij = []
    d_ij = []
    j_atom = []
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


def extract_matrices(file_path):
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
    invlatt = np.linalg.inv(lattice)
    s_ij = []
    for idx, vector in enumerate(r_ij):
        r_i = np.dot(lattice, basis[i_atom-1])
        r_j = np.dot(lattice, basis[j_atoms[idx]-1])
        r_shift = vector/alat
        v_ij = r_shift - r_j - r_i 
        s_ij.append(np.dot(invlatt, v_ij))
        
    return np.array(s_ij, dtype=np.float32).round(4)

# Example usage
file_path = "spin-001/out-3"

alat, lattice = extract_bravais_lattice_matrix(file_path)
matrices = extract_matrices(file_path)
i_atom, r_i, j_atoms, d_ij, r_ij = extract_vectors(file_path)
basis = extract_basis_vectors(file_path)
invbrav = np.linalg.inv(lattice)

print("Lattice:")
print(lattice)
print("Basis:")
print(basis)

end_idx = -1
s_ij = convert_to_maptype_three(i_atom, j_atoms[:end_idx], basis, lattice, alat, r_ij[:end_idx])

outmap = []
for idx, vector in enumerate(r_ij[:end_idx]):
    i_list = [i_atom, j_atoms[idx], s_ij[idx], matrices[idx].flatten(), d_ij[idx]/alat]
    outmap.append(flatten_and_concatenate(i_list))
outmap = np.array(outmap)

fmt = "%d %d   %f %f %f   %f %f %f  %f %f %f %f %f %f     %f"
np.savetxt("outmap.txt", outmap, fmt=fmt)
for idx, vector in enumerate(r_ij[:end_idx]):
    print("|----> ", idx, vector/alat)
    print(" |---> ", basis[i_atom-1])
    print("  |-->", basis[j_atoms[idx]-1])
    #r_ij_mt = np.dot(invbrav, vector/alat).round(2)
    print("   |->",i_atom, j_atoms[idx], ':', s_ij[idx], "|", (matrices[idx]).flatten().round(3))
    #print(i_atom, j_atoms[idx], d_ij[idx], vector, "|", (matrices[idx]).flatten())
    # print(i_atom, j_atoms[idx], d_ij[idx], vector, '|', transform_matrix(matrices[idx]).flatten())
