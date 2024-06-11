#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <Python.h>


static PyObject* process_lines(PyObject* self, PyObject* args) {
    import_array();

    PyObject* input;
    int n_atoms;

    if (!PyArg_ParseTuple(args, "Oi", &input, &n_atoms)) {
        return NULL;  // Error parsing arguments
    }

    PyObject* lines = PyList_GetSlice(input, 0, n_atoms);
    if (lines == NULL) {
        return NULL;  // Error getting list slice
    }

    npy_intp dims[2] = {n_atoms, 3};
    PyObject* xyz = PyArray_SimpleNew(2, dims, NPY_FLOAT);
    if (xyz == 0) {
        PyErr_SetString(PyExc_ValueError, "Could not create numpy array");
        return NULL;
    }

    PyObject* atoms = PyList_New(n_atoms);
    if (atoms == NULL) {
        PyErr_SetString(PyExc_ValueError, "Could not create atoms list");
        return NULL;
    }

    for (int i = 0; i < n_atoms; i++) {
        PyObject* line = PyList_GetItem(lines, i);
        if (line == NULL) {
            PyErr_SetString(PyExc_ValueError, "Could not get line from list");
            return NULL;
        }

        char* line_str = PyUnicode_AsUTF8(line);
        if (line_str == NULL) {
            PyErr_SetString(PyExc_ValueError, "Could not convert line to string");
            return NULL;
        }

        char atom[5];
        float x, y, z;
        int ret = sscanf(line_str, "%s %f %f %f", atom, &x, &y, &z);
        if (ret != 4) {
            PyErr_SetString(PyExc_ValueError, "Could not parse line");
            return NULL;
        }

        float* xyz_line = (float*)PyArray_GETPTR2((PyArrayObject*)xyz, i, 0);
        if (xyz_line == NULL) {
            PyErr_SetString(PyExc_ValueError, "Could not get pointer to numpy array row");
            return NULL;
        }
        xyz_line[0] = x;
        xyz_line[1] = y;
        xyz_line[2] = z;

        PyObject* atom_name = Py_BuildValue("s", atom);
        if (atom_name == NULL) {
            PyErr_SetString(PyExc_ValueError, "Could not create atom name string");
            return NULL;
        }
        PyList_SetItem(atoms, i, atom_name);
    }

    Py_DECREF(lines);

    return Py_BuildValue("[O,O]", atoms, xyz);
}

static PyMethodDef MyTestMethods[] = {
    {"process_lines", process_lines, METH_VARARGS, "Process lines from frame string."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mytestmodule = {
    PyModuleDef_HEAD_INIT,
    "mytest",
    NULL,
    -1,
    MyTestMethods
};

PyMODINIT_FUNC PyInit_mytest(void) {
    return PyModule_Create(&mytestmodule);
}