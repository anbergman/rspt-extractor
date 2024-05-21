"""
This module contains a class, RsptExchange, which is used for extracting and processing
exchange-related data from RSPt calculations.
The class provides methods for extracting exchange data from an input file,
processing the extracted data, and saving the processed data to an output file.
"""

import numpy as np
import rspt_extract as rs


class RsptExchange:
    """
    Class for extracting and processing exchange data from a file.


    The RsptExchange class has the following attributes:
    - file_path: Path to the input file.
    - alat: Lattice constant.
    - lattice: Bravais lattice matrix.
    - matrices: Exchange matrices.
    - i_atom: Atom indices.
    - r_i: Position vectors of the i-th atom.
    - j_atoms: Atom indices of neighboring atoms.
    - d_ij: Distance vectors between atoms.
    - r_ij: Distance vectors between atoms in reduced coordinates.
    - basis: Basis vectors.
    - s_ij: Maptype three representation of distance vectors.
    - outmap: Processed data.

    The RsptExchange class has the following methods:
    - extract_data(): Extracts exchange data from the input file.
    - process_data(): Processes the extracted data.
    - save_output(output_path): Saves the processed data to a file.

    Example usage:
    file_path = "spin-001/out-3"
    output_path = "outmap.txt"

    rspt_exchange = RsptExchange(file_path)
    rspt_exchange.extract_data()
    rspt_exchange.process_data()
    rspt_exchange.save_output(output_path)
    """

    def __init__(self, file_path):
        """
        Initialize the RsptExchange object.

        Args:
            file_path (str): Path to the input file.
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
        Extract exchange data from the input file.
        """
        self.alat, self.lattice = rs.extract_bravais_lattice_matrix(self.file_path)
        self.matrices = rs.extract_exchange_matrices(self.file_path)
        self.i_atom, self.r_i, self.j_atoms, self.d_ij, self.r_ij = (
            rs.extract_distance_vectors(self.file_path)
        )
        self.basis = rs.extract_basis_vectors(self.file_path)
        self.s_ij = rs.convert_to_maptype_three(
            self.i_atom, self.j_atoms, self.basis, self.lattice, self.alat, self.r_ij
        )

    def process_data(self):
        """
        Process the extracted data.
        """
        outmap = []
        for idx, _ in enumerate(self.r_ij):
            i_list = [
                self.i_atom,
                self.j_atoms[idx],
                self.s_ij[idx],
                rs.transform_matrix(self.matrices[idx]).flatten(),
                self.d_ij[idx] / self.alat,
            ]
            outmap.append(rs.flatten_and_concatenate(i_list))
        self.outmap = np.array(outmap)

    def save_output(self, output_path):
        """
        Save the processed data to a file.

        Args:
            output_path (str): Path to the output file.
        """
        fmt1 = "%4d %4d   % 4.1f % 4.1f % 4.1f   "
        fmt2 = "% 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f  % 10.6f % 10.6f % 10.6f     %8.4f"
        fmt = fmt1 + fmt2
        np.savetxt(output_path, self.outmap, fmt=fmt)

def extract_projections(exchange):
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
