#!/usr/bin/env python
# coding: utf-8
'''
ham_check.py
Reads UppASD structure file and reports on missing interactions
and problematic mappings
'''

# Import numpy
import numpy as np

class BColors:
    '''
    Class containing ANSI color codes for colored printing
    Credit: https://stackoverflow.com/questions/287871
    '''
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

############################################################
# Read atom positions from UppASD position file
############################################################
def read_posfile(posfile):
    '''
    Read atom positions from UppASD position file
    '''
    with open(posfile, 'r', encoding='utf-8') as pfile:
        lines = pfile.readlines()
        _positions = np.empty([0, 3])
        _numbers = []
        for line in lines:
            line_data = line.rstrip('\n').split()
            if len(line_data) > 0:
                _positions = np.vstack(
                    (_positions, np.asarray(line_data[2:5]).astype(np.float64)))
                _numbers = np.append(_numbers, np.asarray(
                    line_data[1]).astype(np.int32))

    return _positions, _numbers


############################################################
# Read important keywords from UppASD inputfile `inpsd.dat`
############################################################
def read_inpsd(ifile):
    '''
    Read important keywords from UppASD inputfile `inpsd.dat`
    '''
    _posfiletype = 'C'
    with open(ifile, 'r', encoding='utf-8') as _infile:
        lines = _infile.readlines()
        for idx, line in enumerate(lines):
            line_data = line.rstrip('\n').split()
            if len(line_data) > 0:
                # Find the simulation id
                if line_data[0] == 'simid':
                    _simid = line_data[1]

                # Find the cell data
                if line_data[0] == 'cell':
                    cell = []
                    _lattice = np.empty([0, 3])
                    line_data = lines[idx+0].split()
                    cell = np.append(cell, np.asarray(line_data[1:4]))
                    _lattice = np.vstack(
                        (_lattice, np.asarray(line_data[1:4])))
                    line_data = lines[idx+1].split()
                    cell = np.append(cell, np.asarray(line_data[0:3]))
                    _lattice = np.vstack(
                        (_lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx+2].split()
                    cell = np.append(cell, np.asarray(line_data[0:3]))
                    _lattice = np.vstack(
                        (_lattice, np.asarray(line_data[0:3])))
                    _lattice = _lattice.astype(np.float64)

                # Find the size of the simulated cell
                if line_data[0] == 'ncell':
                    ncell_x = int(line_data[1])
                    ncell_y = int(line_data[2])
                    ncell_z = int(line_data[3])
                    _mesh = [ncell_x, ncell_y, ncell_z]

                # Read the name of the position file
                if line_data[0].strip() == 'posfile':
                    _positions, _numbers = read_posfile(line_data[1])

                # Read the type of coordinate representation
                if line_data[0].strip() == 'posfiletype':
                    _posfiletype = line_data[1]

    return _lattice, _positions, _numbers, _simid, _mesh, _posfiletype


def print_found(_ia, _ja, jij, jji, checkstr):
    '''Print message if coupling found'''
    if checkstr == ' Check ok!':
        s_err = checkstr
    else:
        s_err = BColors.FAIL + checkstr + BColors.ENDC

    print(
        f'ia: {_ia:6.0f} ja: {_ja:6.0f}  Jij: {jij:10.5f}    Jji: {jji:10.5f}  ' + s_err
    )


def print_not_found(_ia, _ja, jij, checkstr):
    '''Print message if coupling not found'''
    s_err = BColors.WARNING + '   Jji:   Not found!' + BColors.FAIL + f'{checkstr:s}' + BColors.ENDC
    print( f'ia: {_ia:6.0f} ja: {_ja:6.0f}  Jij: {jij:10.5f}',s_err)


def check_couplings(ham, n_atoms):
    ''' Check UppASD Hamiltonian for missing couplings'''
    # Reset counters
    i_ok = 0
    i_err = 0

    # Loop over the six first atoms
    for i in range(1, n_atoms+1):
        ham_i = ham[np.int16(ham[:, 0]) == i]
        print('Checking Hamiltonian mapping for atom ', i,
              '    (', ham_i.shape[0], ' interactions)')
        # For every atom, check if Jij = Jji or not
        for row in ham_i:
            _ia = np.int16(row[0])
            _ja = np.int16(row[1])
            jij = row[7]
            ham_j = ham[np.int16(ham[:, 0]) == _ja]
            jji = (ham_j[ham_j[:, 1] == _ia, 7])
            if jji.shape[0] == 1:
                # print('ia:',ia,' ja:',ja,'    Jij:',jij,'    Jji:',jji)
                if jij == jji[0]:
                    checkstr = ' Check ok!'
                    i_ok = i_ok+1
                else:
                    checkstr = ' Check failed!'
                    i_err = i_err+1
                # print('ia: %4d ja: %4d    Jij: % 10.5f    Jji: % 10.5f    %s' %
                #      (_ia, _ja, jij, jji, checkstr))
                print_found(_ia, _ja, jij, jji[0], checkstr)
                # print('ia: %8d ja: %8d    Jij: % 10.5f    Jji: % 10.5f' %(ja,ia,jji,jij))
            else:
                checkstr = ' Check failed!'
                i_err = i_err+1
                print_not_found(_ia, _ja, jij,  checkstr)
                # print('ia: %8d ja: %8d    Jij: % 10.5f    Jji  not found!    %s' %
                #      (_ia, _ja, jij, checkstr))

    if i_err>0:
        s_err_c = BColors.FAIL + f' incorrect mappings: {i_err:5.0f}' + BColors.ENDC
    else:
        s_err_c = BColors.OKGREEN + f' incorrect mappings: {i_err:5.0f}' + BColors.ENDC

    print('Total number of interactions:', i_err+i_ok,
          ', correct mappings: ', i_ok, s_err_c)


############################################################
# Open and read input files
############################################################
INFILE = 'inpsd.dat'
lattice, positions, numbers, simid, mesh, posfiletype = read_inpsd(INFILE)

# Load struct file (Hamiltonian)
hamiltonian = np.genfromtxt('struct.'+simid+'.out')

# Check couplings for atoms in unit cell
check_couplings(hamiltonian, len(positions))
