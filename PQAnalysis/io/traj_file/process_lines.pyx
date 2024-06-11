"""
This module contains the functions that are used to parse the body of 
xyz formatted files. The functions are written in Cython and are compiled
to C code.
"""

import numpy as np
cimport numpy as np
from libc.stdio cimport sscanf

def process_lines_with_atoms(input, int n_atoms):
    """
    Function to parse the body of an xyz file. The function expects a list of
    strings, where each string is a line in the file. The function will parse
    the lines and return a tuple with two elements. The first element is a list
    of strings, where each string is the atom type of the corresponding atom in
    the xyz file. The second element is a numpy array with the xyz coordinates
    of the atoms in the file.

    Parameters
    ----------
    input : list
        A list of strings, where each string is a line in the xyz file.
    n_atoms : int
        The number of atoms in the xyz file.

    Returns
    -------
    tuple
        A tuple with two elements. The first element is a list of strings, where
        each string is the atom type of the corresponding atom in the xyz file.
        The second element is a numpy array with the xyz coordinates of the atoms
        in the file.
    """
    
    cdef np.ndarray[np.float32_t, ndim=2] xyz = np.zeros((n_atoms, 3), dtype=np.float32)
    cdef list atoms = [None] * n_atoms
    cdef char[5] atom
    cdef float x, y, z
    cdef int ret

    for i in range(n_atoms):
        line = input[i]
        line_str = line.encode('utf-8')

        ret = sscanf(line_str, "%s %f %f %f", atom, &x, &y, &z)
        if ret != 4:
            raise ValueError("Could not parse line")

        xyz[i, 0] = x
        xyz[i, 1] = y
        xyz[i, 2] = z

        atoms[i] = atom.decode('utf-8')

    return atoms, xyz

def process_lines(input, int n_atoms):
    """
    Function to parse the body of an xyz file. The function expects a list of
    strings, where each string is a line in the file. The function will parse
    the lines and return a numpy array with the xyz coordinates of the atoms
    in the file.

    Parameters
    ----------
    input : list
        A list of strings, where each string is a line in the xyz file.
    n_atoms : int
        The number of atoms in the xyz file.

    Returns
    -------
    numpy.ndarray
        A numpy array with the xyz coordinates of the atoms in the file.
    """

    cdef np.ndarray[np.float32_t, ndim=2] xyz = np.zeros((n_atoms, 3), dtype=np.float32)
    cdef float x, y, z
    cdef int ret

    for i in range(n_atoms):
        line = input[i]
        line_str = line.encode('utf-8')

        ret = sscanf(line_str, "%*s %f %f %f", &x, &y, &z)
        if ret != 3:
            raise ValueError("Could not parse line")

        xyz[i, 0] = x
        xyz[i, 1] = y
        xyz[i, 2] = z

    return xyz