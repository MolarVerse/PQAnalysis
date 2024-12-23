#ifndef PROCESS_LINES_HPP
#define PROCESS_LINES_HPP

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

namespace py = pybind11;

std::pair<std::vector<std::string>, py::array_t<float>> process_lines_with_atoms(
    const std::vector<std::string>& input,
    int                             n_atoms
);

#endif