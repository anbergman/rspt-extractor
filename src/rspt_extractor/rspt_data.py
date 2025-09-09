"""
This module contains a class, RsptData, which is used for extracting and processing
data from RSPt calculations.
The class provides methods for extracting exchange data from an input file,
processing the extracted data, and saving the processed data to an output file.
"""

import importlib.resources
import numpy as np
import spglib as spg
from .rspt_extract import (
    check_relativistic,
    extract_moments_scf,
    extract_position_scf,
    extract_position_data,
    extract_lattice_scf,
)

# Simple periodic table: index = atomic number, value = element symbol
ELEMENTS = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]


class RsptData:
    """
    A class to handle RSPT SCF data extraction and file generation.

    This class provides methods to extract data from an input file and generate
    various output files required for UppASD simulations. It handles lattice
    parameters, basis vectors, and magnetic moments.

    Attributes:
        alat (float): Lattice constant.
        lattice (np.ndarray): Bravais lattice matrix.
        basis (np.ndarray): Basis vectors.
        moments (np.ndarray): Magnetic moments.

    Methods:
        extract_scf_data():
        print_lattice(output_path):
            Print the lattice matrix to a text file.
        print_positions(output_path):
            Print the basis vectors to an UppASD posfile.
        print_moments(output_path):
            Print the magnetic moments to an UppASD posfile.
        print_template(posfile, momfile, exchange, file_path, mtype):

    Author: Anders Bergman
    """

    def __init__(self, args):
        """
        Initializes the class with the given file path and extracts data.

        Args:
            file_path (str): The path to the file containing the data.

        Attributes:
            file_path (str): The path to the file containing the data.
            alat (float or None): The lattice constant, initialized as None.
            lattice (list or None): The lattice vectors, initialized as None.
            basis (list or None): The basis atoms, initialized as None.
            moments (list or None): The magnetic moments, initialized as None.

        Methods:
            extract_scf_data(): Extracts data from the file and populates the
                            attributes.

        Author: Anders Bergman
        """
        self.scf_file = args.scf
        self.data_file = args.data
        self.alat = None
        self.lattice = None
        self.basis = None
        self.moments = None
        self.atoms = None
        self.types = []
        self.species = []
        self.filter_list = []
        self.is_relativistic = False
        self.threshold = args.threshold
        self.bohr_to_angstrom = 0.529177249
        self.extract_scf_data()
        self.get_symmetry_data()

    def extract_scf_data(self):
        """
        Extracts and sets the lattice matrix, basis vectors, and moments from
        the specified file.

        This method reads the file located at `self.file_name` and extracts
        the Bravais lattice matrix, basis vectors, and magnetic moments using
        the `rs` module's extraction functions. The extracted data is then
        assigned to the instance variables `self.alat`, `self.lattice`,
        `self.basis`, and `self.moments`.

        Args:
            None

        Returns:
            None

        Raises:
            IOError: If the file at `self.file_name` cannot be read.
            ValueError: If the extracted data is invalid or incomplete.

        Author:
            Anders Bergman
        """
        # print(50 * "-")
        scf_moments = []
        species = []
        # Check if the SCF file is provided and if the
        # calculation is relativistic or not
        # Extract lattice and moments from the SCF file
        if self.scf_file:
            self.is_relativistic = check_relativistic(self.scf_file)
            # print(f"Is the calculation relativistic? {self.is_relativistic}")

            self.alat, self.lattice = extract_lattice_scf(self.scf_file)
            print("Input lattice matrix: [alat]")
            for row in self.lattice:
                print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
            print(f"Input lattice constant [a.u.]:\n  {self.alat}")
            alat_angstrom = self.alat * 0.529177
            print(f"Input lattice constant [Å]:\n  {alat_angstrom}")
            scf_moments = extract_moments_scf(self.scf_file, self.is_relativistic)

        # Extract basis vectors from the data file
        if self.scf_file:
            self.types, self.basis, species = extract_position_scf(self.scf_file)
            self.atoms = range(0, len(self.types))
            # Ensure all basis vectors are within 0 < x < 1
            # self.basis = np.mod(self.basis, 1.0)
            print("Input basis vectors: (crystal)")
            for row in self.basis:
                print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")

        # # Extract basis vectors from the data file
        elif self.data_file:
            self.types, self.basis, species = extract_position_data(self.data_file)
            self.atoms = range(0, len(self.types))
            # Ensure all basis vectors are within 0 < x < 1
            # self.basis = np.mod(self.basis, 1.0)
            print("Input basis vectors: (crystal)")
            for row in self.basis:
                print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
        # Extract moments from the SCF file
        # if self.scf_file and self.data_file:
        if self.scf_file:
            self.moments = []
            self.species = []
            for atom in self.types:
                # print(f"Atom type: {atom}")
                self.moments.append(scf_moments[atom - 1].tolist())
                self.species.append(species[atom - 1])

            # Print species
            print(f"Input species list: {self.species}")

            # Print moments
            total_moments = np.array([mom[2] for mom in self.moments])
            print("SCF magnetic moments: [mu_B/atom] ")
            print(total_moments)

            # Filter out moments with magnitude below the threshold
            print(
                f"Moment threshold: {self.threshold} : {total_moments>=self.threshold}"
            )
            tmp_moments = []
            tmp_basis = []
            tmp_atoms = []
            for i, moment in enumerate(self.moments):
                # Check total (mS+mL) moment
                if np.abs(moment[2]) >= self.threshold:
                    tmp_moments.append(moment)
                    tmp_basis.append(self.basis[i])
                    tmp_atoms.append(self.atoms[i])
                    self.filter_list.append(i + 1)

            self.full_moments = np.array(self.moments)
            self.moments = np.array(tmp_moments)
            self.full_basis = np.array(self.basis)
            self.basis = np.array(tmp_basis)
            self.full_atoms = np.array(self.atoms)
            self.atoms = np.array(tmp_atoms)
            del tmp_moments
            del tmp_basis
            del tmp_atoms

    def print_lattice(self, output_path):
        """
        Save the lattice configuration to a file.

        This method saves the lattice configuration to a specified file path
        using a fixed format for floating-point numbers.

        Args:
            output_path (str): The path to the file where the lattice
                       configuration will be saved.

        Returns:
            None

        Example:
            >>> obj.print_lattice('/path/to/output.txt')

        Author: Anders Bergman
        """
        np.savetxt(output_path, self.lattice, fmt="% 10.6f % 10.6f % 10.6f")

    def print_positions(self, output_path):
        """
        Save atomic positions to a file.

        This method iterates over the atomic basis positions, formats them, and
        saves them to a specified output file.

        Args:
            output_path (str): The path to the output file where the positions
                will be saved.

        Returns:
            None

        Example:
            >>> obj.print_positions('output.txt')

        Author:
            Anders Bergman
        """
        posdata = []
        for i, pos in enumerate(self.basis):
            posdata.append(np.hstack([i + 1, i + 1, pos]))
        np.savetxt(
            output_path, np.array(posdata), fmt="%4d %4d  % 10.6f % 10.6f % 10.6f"
        )

    def print_poscar(self, output_path, comment=None):
        """
        Write a VASP POSCAR file using direct (fractional) coordinates.

        The lattice vectors are converted from atomic units (Bohr) to
        Angstrom using `self.bohr_to_angstrom` and written as the POSCAR
        lattice. The basis positions in `self.basis` are assumed to be
        fractional (direct) coordinates and are written under the
        `Direct` label. Species and counts are inferred from the instance's
        species/atoms data when available; otherwise a generic element
        label `X` is used.

        Args:
            output_path (str): Path to the output POSCAR file.
            comment (str, optional): Optional comment/label for the POSCAR
                first line. If not provided a default header is written.

        Returns:
            None
        """
        # Prepare header
        header = (
            comment
            if comment is not None
            else f"RsptData POSCAR (alat={self.alat * self.bohr_to_angstrom:.3f} Å)"
        )

        # Prepare atom species for the (possibly filtered) atom list
        # Determine which atom indices are currently active
        if hasattr(self, "atoms") and self.atoms is not None:
            try:
                atom_indices = list(map(int, list(self.atoms)))
            except Exception:
                atom_indices = list(range(len(self.basis)))
        else:
            atom_indices = list(range(len(self.basis)))

        species_for_atoms = None
        if getattr(self, "species", None):
            try:
                species_for_atoms = [self.species[i] for i in atom_indices]
            except Exception:
                species_for_atoms = None

        # If species are not available, fallback to single generic element
        if not species_for_atoms:
            nat = len(self.basis)
            element_list = ["X"]
            counts = [nat]
            species_for_atoms = ["X"] * nat
        else:
            # Convert species identifiers to element symbols and count them
            mapped = [self._species_to_symbol(s) for s in species_for_atoms]
            element_list = []
            counts = []
            for s in mapped:
                if s not in element_list:
                    element_list.append(s)
                    counts.append(1)
                else:
                    counts[element_list.index(s)] += 1

        # Build POSCAR text lines
        lines = []
        lines.append(header)
        lines.append(str(self.alat * self.bohr_to_angstrom))
        for vec in self.lattice:
            lines.append(f"{vec[0]:16.10f} {vec[1]:16.10f} {vec[2]:16.10f}")
        lines.append(" ".join(element_list))
        lines.append(" ".join(str(c) for c in counts))
        lines.append("Direct")

        # Basis positions are assumed to be direct/fractional
        for pos in np.array(self.basis, dtype=float):
            lines.append(f"{pos[0]:16.10f} {pos[1]:16.10f} {pos[2]:16.10f}")

        # Write file
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")

    def _species_to_symbol(self, s):
        """
        Convert a species identifier to an element symbol.

        Accepts integer atomic numbers (1..Z) or string labels and returns a
        chemical symbol. Falls back to the string representation if unknown.
        """
        # If it's an integer or can be interpreted as one, map via ELEMENTS
        try:
            ani = int(s)
            if 0 <= ani < len(ELEMENTS):
                return ELEMENTS[ani]
        except Exception:
            pass

        # If it's already a string like 'Fe' or '26', return cleaned string
        s_str = str(s).strip()
        # If it looks like a number, try mapping
        if s_str.isdigit():
            idx = int(s_str)
            if 0 <= idx < len(ELEMENTS):
                return ELEMENTS[idx]

        # Otherwise return as-is
        return s_str

    def print_moments(self, output_path):
        """
        Saves the magnetic moments to a specified file.

        This method iterates over the magnetic moments stored in the `moments`
        attribute, formats them, and saves them to a text file at the given
        `output_path`.

        Args:
            output_path (str): The path to the file where the moments will be saved.

        Returns:
            None

        Example:
            >>> obj.print_moments('output.txt')

        Notes:
            The output file will contain columns for the index, a constant value
            of 1, the x-component of the moment, and three additional constant
            values (0.0, 0.0, 1.0).

        Author: Anders Bergman
        """
        momdata = []
        for i, mom in enumerate(self.moments):
            momdata.append(np.hstack([i + 1, 1, mom[0], 0.0, 0.0, 1.0]))
        np.savetxt(
            output_path,
            np.array(momdata),
            fmt="%4d %4d   % 10.6f   % 10.6f % 10.6f % 10.6f",
        )

    def print_template(self, posfile, momfile, exchange, file_path, mtype):
        """
        Generates and writes a formatted output file based on a template and
        provided parameters.

        This method reads a template file, substitutes placeholders with
        corresponding values from the provided parameters and the object's
        lattice attribute, and writes the formatted output to a specified file.

        Args:
            posfile (str): Path to the position file.
            momfile (str): Path to the moment file.
            exchange (str): Exchange parameter value.
            file_path (str): Path where the output file will be saved.

        Raises:
            FileNotFoundError: If the template file is not found.
            IOError: If there is an error reading the template file or writing
                     the output file.

        Author: Anders Bergman
        """
        # Digest proper maptype and posfiletype from mtype
        if mtype == "3":
            maptype = "3"
            posfiletype = "D"
        elif mtype == "D":
            maptype = "1"
            posfiletype = "D"
        elif mtype == "C":
            maptype = "1"
            posfiletype = "C"
        # Create a dictionary with all placeholders and their values
        template_values = {
            "cell_00": self.lattice[0][0],
            "cell_01": self.lattice[0][1],
            "cell_02": self.lattice[0][2],
            "cell_10": self.lattice[1][0],
            "cell_11": self.lattice[1][1],
            "cell_12": self.lattice[1][2],
            "cell_20": self.lattice[2][0],
            "cell_21": self.lattice[2][1],
            "cell_22": self.lattice[2][2],
            "posfile": posfile,
            "momfile": momfile,
            "exchange": exchange,
            "maptype": maptype,
            "posfiletype": posfiletype,
            "alat": self.alat * self.bohr_to_angstrom * 1.0e-10,
        }

        # Read the template from the package using importlib.resources
        with importlib.resources.files("rspt_extractor").joinpath(
            "inpsd.template"
        ).open("r", encoding="utf-8") as file:
            template = file.read()
        # Read the template from the external file
        #   with open("inpsd.template", "r", encoding="utf-8") as file:
        #   template = file.read()

        # Use the template to create the final output string
        output = template.format(**template_values)

        # Save the output to the file
        with open(file_path, "w", encoding="utf-8") as file:
            file.write(output)

    def get_symmetry_data(self):
        """
        Get the symmetry data of the system.

        This method uses the spglib library to calculate the symmetry data
        of the system.

        Args:
            None

        Returns:
            None
        """
        print(50 * "-")
        # cell = (self.lattice, self.basis, self.species)
        print("Lattice matrix:")
        for row in self.lattice:
            print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
        print("Basis vectors:")
        for row in self.full_basis:
            print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
        print("Species: ", self.species)
        print("Magnetic moments:")
        print(self.full_moments)

        cell = (self.lattice, self.full_basis, self.species, self.full_moments)
        # Get the symmetry data
        symmetry_data = spg.get_symmetry_dataset(cell=cell)
        # Print the symmetry data
        print("Space group number: ", symmetry_data.number)
        print("Space group symbol: ", symmetry_data.international)
        print("Hall group: ", symmetry_data.hall)
        print("Number of symmetry operations: ", len(symmetry_data.rotations))
        print(50 * "-")


if __name__ == "__main__":
    # Example usage
    FILE_PATH = "spin-001/out-scf"
    rspt_exchange = RsptData(FILE_PATH)
    rspt_exchange.print_lattice("lattice.dat")
    rspt_exchange.print_positions("posfile")
    rspt_exchange.print_moments("momfile")
    rspt_exchange.print_template("posfile", "momfile", "jfile", "inpsd.dat", mtype="3")
