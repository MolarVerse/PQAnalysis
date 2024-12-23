#include <format>

#include "cell.hpp"

PYBIND11_MODULE(cell, m)
{
    m.doc() = "Cell module";
    py::class_<Cell>(m, "Cell")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def(py::init<double, double, double, double, double, double>())
        .def("bouding_edges", &Cell::bouding_edges)
        .def("volume", &Cell::volume)
        .def("is_vacuum", &Cell::is_vacuum)
        .def("image", &Cell::image, py::arg("pos"))
        .def(
            "init_from_box_matrix",
            &Cell::init_from_box_matrix,
            py::arg("box_matrix")
        )
        .def_property(
            "box_matrix",
            &Cell::get_box_matrix,
            &Cell::set_box_matrix
        )
        .def_property(
            "box_lengths",
            &Cell::get_box_lengths,
            &Cell::set_box_lengths
        )
        .def_property(
            "box_angles",
            &Cell::get_box_angles,
            &Cell::set_box_angles
        )
        .def_property("x", &Cell::get_x, &Cell::set_x)
        .def_property("y", &Cell::get_y, &Cell::set_y)
        .def_property("z", &Cell::get_z, &Cell::set_z)
        .def_property("alpha", &Cell::get_alpha, &Cell::set_alpha)
        .def_property("beta", &Cell::get_beta, &Cell::set_beta)
        .def_property("gamma", &Cell::get_gamma, &Cell::set_gamma)
        .def("__eq__", &Cell::operator==)
        .def(
            "__repr__",
            [](Cell &a)
            {
                return std::format(
                    "Cell(x={}, y={}, z={}, alpha={}, beta={}, gamma={})",
                    a.get_x(),
                    a.get_y(),
                    a.get_z(),
                    a.get_alpha(),
                    a.get_beta(),
                    a.get_gamma()
                );
            }
        );
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
