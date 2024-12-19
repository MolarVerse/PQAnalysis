#ifndef CELL_HPP
#define CELL_HPP

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <math.h>

namespace py = pybind11;
using array_d = py::array_t<double>;

Eigen::Map<Eigen::Matrix> numpy_to_eigen(array_d &array);

class CoreCell
{
private:
    double _x, _y, _z, _alpha, _beta, _gamma;
    Eigen::Matrix<double, 3, 3> _box_matrix;
    array_d _box_lengths;
    array_d _box_angles;

public:
    CoreCell();
    CoreCell(double x, double y, double z, double alpha, double beta, double gamma);
    ~CoreCell() = default;

private:
    Matrix<double, 3, 3> _init_box_matrix();

public:
    // 2D array
    Eigen::Matrix bouding_edges();
    double volume();
    bool is_vacuum();
    array_d image(array_d pos);
    CoreCell &init_from_box_matrix(array_d box_matrix);

    // Getters
    array_d get_box_matrix() const
    {
        return _box_matrix;
    }
    array_d get_box_lengths() const
    {
        return _box_lengths;
    }
    array_d get_box_angles() const
    {
        return _box_angles;
    }
    double get_x() const { return _x; }
    double get_y() const { return _y; }
    double get_z() const { return _z; }
    double get_alpha() const { return _alpha; }
    double get_beta() const { return _beta; }
    double get_gamma() const { return _gamma; }

    // Setters
    void set_box_matrix(Eigen::Matrix box_matrix)
    {
        _box_matrix = numpy_to_eigen(box_matrix);
    }
    void set_box_lengths(array_d box_lengths)
    {
        _box_lengths = box_lengths;
    }
    void set_box_angles(array_d box_angles)
    {
        _box_angles = box_angles;
    }
    void set_x(double x) { _x = x; }
    void set_y(double y) { _y = y; }
    void set_z(double z) { _z = z; }
    void set_alpha(double alpha) { _alpha = alpha; }
    void set_beta(double beta) { _beta = beta; }
    void set_gamma(double gamma) { _gamma = gamma; }
};

// Convert NumPy array to Eigen Map
Eigen::Map<Eigen::Matrix> numpy_to_eigen(array_d array)
{
    py::buffer_info info = array.request();

    if (info.ndim != 2)
    {
        throw std::runtime_error("Input array must be 2D");
    }

    return Eigen::Map<Eigen::Matrix>(
        static_cast<double *>(info.ptr),
        info.shape[0], // rows
        info.shape[1]  // cols
    );
}

#endif // CELL_HPP