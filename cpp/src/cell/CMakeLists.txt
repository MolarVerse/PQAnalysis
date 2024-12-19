# Create a library called "cell" which includes the source file "cell.cpp".

# collect all source files
file(GLOB CELL "*.cpp")

# include the header files
include_directories(${PROJECT_SOURCE_DIR}/cpp/include/cell)

# create the library
pybind11_add_module(core_cell ${CELL})

# link the Eigen library to the pybind11 module
target_link_libraries(core_cell PRIVATE Eigen3::Eigen)

# link the library to the pybind11 module
target_compile_definitions(core_cell
                           PRIVATE VERSION_INFO=${PROJECT_VERSION})