#ifndef CELL_HPP
#define CELL_HPP

#define STRINGIFY(x)       #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <math.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

namespace py  = pybind11;
using array_d = py::array_t<double>;

// Convert 2D vector to numpy array
array_d vector_2d_to_array(const std::vector<std::vector<double>> &vec);
std::vector<std::vector<double>> array_to_vector_2d(array_d arr);
std::vector<double>              array_to_vector_1d(array_d arr);

class Cell
{
   private:
    std::vector<double>              _box_lengths, _box_angles;
    std::vector<std::vector<double>> _box_matrix;

   public:
    Cell();
    Cell(double x, double y, double z);
    Cell(double x, double y, double z, double alpha, double beta, double gamma);
    ~Cell() = default;

   private:
    std::vector<std::vector<double>> _setup_box_matrix();

   public:
    std::vector<std::vector<double>> bouding_edges();
    double                           volume();
    bool                             is_vacuum();
    array_d                          image(array_d pos);
    Cell                            &init_from_box_matrix(array_d box_matrix);

    // Operators
    bool operator==(const Cell &other) const
    {
        // Compare box matrix with other box matrix contents
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (_box_matrix[i][j] != other._box_matrix[i][j])
                {
                    return false;
                }
            }
        }
    }

    bool isclose(const Cell &other, double rel_tol = 1e-9, double abs_tol = 0.0)
        const
    {
        // Compare box matrix with other box matrix contents
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (fabs(_box_matrix[i][j] - other._box_matrix[i][j]) >
                    fmax(
                        rel_tol * fmax(
                                      fabs(_box_matrix[i][j]),
                                      fabs(other._box_matrix[i][j])
                                  ),
                        abs_tol
                    ))
                {
                    return false;
                }
            }
        }
        return true;
    }

    // Getters
    array_d get_box_matrix() const
    {
        // Convert 2D vector to numpy array
        return vector_2d_to_array(_box_matrix);
    }
    array_d get_box_lengths() const
    {
        return array_d(_box_lengths.size(), _box_lengths.data());
    }
    array_d get_box_angles() const
    {
        return array_d(_box_angles.size(), _box_angles.data());
    }
    double get_x() const { return _box_lengths[0]; }
    double get_y() const { return _box_lengths[1]; }
    double get_z() const { return _box_lengths[2]; }
    double get_alpha() const { return _box_angles[0]; }
    double get_beta() const { return _box_angles[1]; }
    double get_gamma() const { return _box_angles[2]; }

    // Setters
    void set_box_matrix(array_d box_matrix)
    {
        // Convert numpy array to 2D vector
        _box_matrix = array_to_vector_2d(box_matrix);
    }
    void set_box_lengths(array_d box_lengths)
    {
        // Convert numpy array to vector
        _box_lengths = array_to_vector_1d(box_lengths);
    }
    void set_box_angles(array_d box_angles)
    {
        // Convert numpy array to vector
        _box_angles = array_to_vector_1d(box_angles);
    }
    void set_x(double x) { _box_lengths[0] = x; }
    void set_y(double y) { _box_lengths[1] = y; }
    void set_z(double z) { _box_lengths[2] = z; }
    void set_alpha(double alpha) { _box_angles[0] = alpha; }
    void set_beta(double beta) { _box_angles[1] = beta; }
    void set_gamma(double gamma) { _box_angles[2] = gamma; }
};

#endif   // CELL_HPP