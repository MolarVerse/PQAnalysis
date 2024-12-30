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
            // Cast the input list to a vector of strings
            std::vector<std::string> _input =
                py::cast<std::vector<std::string>>(input);

            // try-catch block to catch any exceptions
            try
            {
                // Process the lines and return the result
                auto [atoms, xyz] = process_lines_with_atoms(_input, n_atoms);

                // Cast the atoms vector to a Python list
                py::array xyz_array = py::cast(xyz);

                // Return the atoms and xyz as a pair
                return std::make_pair(atoms, xyz_array.reshape({n_atoms, 3}));
            }
            catch (const std::exception &e)
            {
                // Raise a Python value error with the exception message
                throw py::value_error(e.what());
            }
        },
        py::arg("input"),
        py::arg("n_atoms")
    );
}