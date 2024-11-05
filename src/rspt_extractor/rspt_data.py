"""
This module contains a class, RsptData, which is used for extracting and processing
data from RSPt calculations.
The class provides methods for extracting exchange data from an input file,
processing the extracted data, and saving the processed data to an output file.
"""

import importlib.resources
import numpy as np
from .rspt_extract import (
    check_relativistic,
    extract_moments_scf,
    extract_position_data,
    extract_lattice_scf,
)


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
        print_template(posfile, momfile, exchange, file_path):

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
        self.filter_list = []
        self.is_relativistic = False
        self.threshold = args.threshold
        self.extract_scf_data()

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
        types = []
        scf_moments = []
        # Check if the SCF file is provided and if the
        # calculation is relativistic or not
        # Extract lattice and moments from the SCF file
        if self.scf_file:
            self.is_relativistic = check_relativistic(self.scf_file)
            # print(f"Is the calculation relativistic? {self.is_relativistic}")

            self.alat, self.lattice = extract_lattice_scf(self.scf_file)
            print("SCF lattice matrix: [alat]")
            for row in self.lattice:
                print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
            print(f"SCF lattice constant [a.u.]:\n  {self.alat}")

            scf_moments = extract_moments_scf(self.scf_file, self.is_relativistic)

        # Extract basis vectors from the data file
        if self.data_file:
            types, self.basis = extract_position_data(self.data_file)
            self.atoms = range(0, len(types))
            print("`data` basis vectors: (crystal)")
            for row in self.basis:
                print(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")

        # Extract moments from the SCF file
        if self.scf_file and self.data_file:
            self.moments = []
            for atom in types:
                # print(f"Atom type: {atom}")
                self.moments.append(scf_moments[atom - 1].tolist())
            print("SCF magnetic moments: [mu_B/atom] ")
            print([mom[2] for mom in self.moments])

            # Filter out moments with magnitude below the threshold
            print(f"Moment threshold: {self.threshold}")
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

            self.moments = np.array(tmp_moments)
            self.basis = np.array(tmp_basis)
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

    def print_template(self, posfile, momfile, exchange, file_path):
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


if __name__ == "__main__":
    # Example usage
    FILE_PATH = "spin-001/out-scf"
    rspt_exchange = RsptData(FILE_PATH)
    rspt_exchange.print_lattice("lattice.dat")
    rspt_exchange.print_positions("posfile")
    rspt_exchange.print_moments("momfile")
    rspt_exchange.print_template("posfile", "momfile", "jfile", "inpsd.minimal")
