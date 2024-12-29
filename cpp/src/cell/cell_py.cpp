#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <format>

#include "cell.hpp"

#define STRINGIFY(x)       #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace py = pybind11;

PYBIND11_MODULE(cell, m)
{
    m.doc() = "Cell module";
    py::class_<Cell>(m, "Cell")
        .def(
            py::init<float, float, float, float, float, float>(),
            py::arg("x") = cbrt(std::numeric_limits<float>::max()) *
                           0.99,   // 0.99 * max float for is_vacuum check
            py::arg("y") = cbrt(std::numeric_limits<float>::max()) *
                           0.99,   // 0.99 * max float for is_vacuum check
            py::arg("z") = cbrt(std::numeric_limits<float>::max()) *
                           0.99,   // 0.99 * max float for is_vacuum check
            py::arg("alpha") = 90.0,
            py::arg("beta")  = 90.0,
            py::arg("gamma") = 90.0
        )
        .def_property_readonly(
            "bounding_edges",
            [](Cell &self)
            {
                py::array edges = py::cast(self.bounding_edges());
                return edges.reshape({-1, 3});
            }
        )
        .def("volume", &Cell::volume)
        .def("is_vacuum", &Cell::is_vacuum)
        .def(
            "image",
            [](Cell &self, py::array_t<float> pos)
            {
                // Convert position to a vector of vectors
                py::buffer_info info           = pos.request();
                auto            original_shape = info.shape;

                // Check if the input is a 1D array
                std::vector<float> pos_vec =
                    py::cast<std::vector<float>>(pos.attr("flatten")());

                // Check if the input is a 2D array
                py::array image_array = py::cast(self.image(pos_vec));
                return image_array.reshape(original_shape);
            },
            py::arg("pos")
        )
        .def(
            "init_from_box_matrix",
            [](Cell &self, py::array_t<float> box_matrix)
            {
                py::buffer_info info = box_matrix.request();
                if (info.ndim != 2)
                {
                    throw std::runtime_error("box_matrix must be a 2D array");
                }
                if (info.shape[0] != 3 || info.shape[1] != 3)
                {
                    throw std::runtime_error("box_matrix must be a 3x3 array");
                }
                std::vector<float> box_matrix_vec(
                    info.shape[0] * info.shape[1]
                );
                return self.init_from_box_matrix(box_matrix_vec);
            },
            py::arg("box_matrix")
        )
        .def_property(
            "box_matrix",
            [](Cell &self)
            {
                py::array box_matrix = py::cast(self.get_box_matrix());
                return box_matrix.reshape({3, 3});
            },
            [](Cell &self, py::array_t<float> box_matrix)
            {
                py::buffer_info info = box_matrix.request();
                if (info.ndim != 2)
                {
                    throw std::runtime_error("box_matrix must be a 2D array");
                }
                if (info.shape[0] != 3 || info.shape[1] != 3)
                {
                    throw std::runtime_error("box_matrix must be a 3x3 array");
                }
                std::vector<float> box_matrix_vec(
                    info.shape[0] * info.shape[1]
                );
                self.set_box_matrix(box_matrix_vec);
            }
        )
        .def_property(
            "box_lengths",
            [](Cell &self)
            {
                py::array box_lengths = py::cast(self.get_box_lengths());
                return box_lengths;
            },
            [](Cell &self, py::array_t<float> box_lengths)
            {
                py::buffer_info info = box_lengths.request();
                if (info.ndim != 1)
                {
                    throw std::runtime_error("box_lengths must be a 1D array");
                }
                if (info.shape[0] != 3)
                {
                    throw std::runtime_error("box_lengths must have 3 elements"
                    );
                }
                std::vector<float> box_lengths_vec(info.shape[0]);
                self.set_box_lengths(box_lengths_vec);
            }
        )
        .def_property(
            "box_angles",
            [](Cell &self)
            {
                py::array box_angles = py::cast(self.get_box_angles());
                return box_angles;
            },
            [](Cell &self, py::array_t<float> box_angles)
            {
                py::buffer_info info = box_angles.request();
                if (info.ndim != 1)
                {
                    throw std::runtime_error("box_angles must be a 1D array");
                }
                if (info.shape[0] != 3)
                {
                    throw std::runtime_error("box_angles must have 3 elements");
                }
                std::vector<float> box_angles_vec(info.shape[0]);
                self.set_box_angles(box_angles_vec);
            }
        )
        .def_property("x", &Cell::get_x, &Cell::set_x)
        .def_property("y", &Cell::get_y, &Cell::set_y)
        .def_property("z", &Cell::get_z, &Cell::set_z)
        .def_property("alpha", &Cell::get_alpha, &Cell::set_alpha)
        .def_property("beta", &Cell::get_beta, &Cell::set_beta)
        .def_property("gamma", &Cell::get_gamma, &Cell::set_gamma)
        .def(
            "isclose",
            &Cell::isclose,
            py::arg("other"),
            py::arg("rtol") = 1e-9,
            py::arg("atol") = 0.0,
            "Check if two cells are close"
        )
        .def("__eq__", [](Cell &a, Cell &b) { return a.isclose(b); })
        .def("__str__", &Cell::to_string)
        .def("__repr__", &Cell::to_string);
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
