#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "process_lines.hpp"

namespace py = pybind11;

PYBIND11_MODULE(process_lines, m)
{
    m.doc() = "Process XYZ file lines using PyBind11";
    m.def(
        "process_lines_with_atoms",
        [](const py::list &input, int n_atoms)
        {
            std::vector<std::string> _input =
                py::cast<std::vector<std::string>>(input);
            auto [atoms, xyz]    = process_lines_with_atoms(_input, n_atoms);
            py::array xyz_array  = py::cast(xyz);
            py::list  atoms_list = py::cast(atoms);

            return std::make_pair(atoms_list, xyz_array.reshape({n_atoms, 3}));
        },
        py::arg("input"),
        py::arg("n_atoms")
    );
}