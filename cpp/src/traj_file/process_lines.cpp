#include "process_lines.hpp"

#include <sstream>

std::pair<std::vector<std::string>, py::array_t<float>> process_lines_with_atoms(
    const std::vector<std::string>& input,
    int                             n_atoms
)
{
    std::vector<std::string> atoms(n_atoms);
    auto                     xyz         = py::array_t<float>({n_atoms, 3});
    auto                     xyz_mutable = xyz.mutable_unchecked<2>();

    for (int i = 0; i < n_atoms; i++)
    {
        std::istringstream iss(input[i]);
        std::string        atom;
        float              x, y, z;

        if (!(iss >> atom >> x >> y >> z))
        {
            throw std::runtime_error("Could not parse line");
        }

        atoms[i]          = atom;
        xyz_mutable(i, 0) = x;
        xyz_mutable(i, 1) = y;
        xyz_mutable(i, 2) = z;
    }

    return {atoms, xyz};
}