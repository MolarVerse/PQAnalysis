#include "cell.hpp"
#include <format>

PYBIND11_MODULE(core_cell, m)
{
    m.doc() = "Cell module";
    py::class_<CoreCell>(m, "CoreCell")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double, double>())
        .def("bouding_edges", &CoreCell::bouding_edges)
        .def("volume", &CoreCell::volume)
        .def("is_vacuum", &CoreCell::is_vacuum)
        .def("image", &CoreCell::image, py::arg("pos"))
        .def("init_from_box_matrix", &CoreCell::init_from_box_matrix, py::arg("box_matrix"))
        .def_property("box_matrix", &CoreCell::get_box_matrix, &CoreCell::set_box_matrix)
        .def_property("box_lengths", &CoreCell::get_box_lengths, &CoreCell::set_box_lengths)
        .def_property("box_angles", &CoreCell::get_box_angles, &CoreCell::set_box_angles)
        .def_property("x", &CoreCell::get_x, &CoreCell::set_x)
        .def_property("y", &CoreCell::get_y, &CoreCell::set_y)
        .def_property("z", &CoreCell::get_z, &CoreCell::set_z)
        .def_property("alpha", &CoreCell::get_alpha, &CoreCell::set_alpha)
        .def_property("beta", &CoreCell::get_beta, &CoreCell::set_beta)
        .def_property("gamma", &CoreCell::get_gamma, &CoreCell::set_gamma)
        .def("__repr__",
             [](CoreCell &a)
             {
                 return std::format("CoreCell(x={}, y={}, z={}, alpha={}, beta={}, gamma={})",
                                    a.get_x(), a.get_y(), a.get_z(), a.get_alpha(), a.get_beta(), a.get_gamma());
             });
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
