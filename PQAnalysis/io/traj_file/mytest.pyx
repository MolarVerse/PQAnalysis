# mytest.pyx
import numpy as np
cimport numpy as np
from libc.stdio cimport sscanf

def process_lines(input, int n_atoms):
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