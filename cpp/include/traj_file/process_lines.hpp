#ifndef PROCESS_LINES_HPP
#define PROCESS_LINES_HPP

#include <string>
#include <vector>

std::pair<std::vector<std::string>, std::vector<float>> process_lines_with_atoms(
    const std::vector<std::string>& input,
    int                n_atoms
);

#endif