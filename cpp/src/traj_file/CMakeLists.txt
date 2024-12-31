# Create a library called "traj_file" which includes the source file "traj_file.cpp".

# collect all source files
file(GLOB PROCESS_LINES "*.cpp")

# include the header files
include_directories(${PROJECT_SOURCE_DIR}/cpp/include/traj_file)

# create the library
pybind11_add_module(process_lines ${PROCESS_LINES})

# link the library to the pybind11 module
target_compile_definitions(process_lines
                           PRIVATE VERSION_INFO=${VERSION_INFO})
# link openmp
target_link_libraries(process_lines PRIVATE OpenMP::OpenMP_CXX)