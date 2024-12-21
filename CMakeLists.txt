cmake_minimum_required(VERSION 3.4...3.18)
project(PQAnalysis)

# set the C++ standard
set(CMAKE_CXX_STANDARD 20)

# find pybind11 package
find_package(pybind11 REQUIRED)

# add cpp subdirectory
add_subdirectory(cpp/src)