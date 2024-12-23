#include "process_lines.hpp"

PYBIND11_MODULE(process_lines, m)
{
    m.doc() = "Process XYZ file lines using PyBind11";
    m.def(
        "process_lines_with_atoms",
        &process_lines_with_atoms,
        "Process xyz file lines",
        py::arg("input"),
        py::arg("n_atoms")
    );
}