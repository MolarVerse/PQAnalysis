#ifndef CELL_HPP
#define CELL_HPP

#include <math.h>

#include <vector>

class Cell
{
   private:
    std::vector<double>              _box_lengths, _box_angles;
    std::vector<std::vector<double>> _box_matrix;

   public:
    Cell();
    Cell(double x, double y, double z, double alpha, double beta, double gamma);
    ~Cell() = default;

   private:
    std::vector<std::vector<double>> _setup_box_matrix();

   public:
    std::vector<std::vector<double>> bounding_edges();
    double                           volume();
    bool                             is_vacuum();
    std::vector<std::vector<double>> image(std::vector<std::vector<double>> pos
    );
    Cell &init_from_box_matrix(std::vector<std::vector<double>> box_matrix);

    bool isclose(const Cell &other, double rtol = 1e-9, double atol = 0.0) const
    {
        // Compare box matrix with other box matrix contents
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (fabs(_box_matrix[i][j] - other._box_matrix[i][j]) >
                    fmax(
                        rtol * fmax(
                                   fabs(_box_matrix[i][j]),
                                   fabs(other._box_matrix[i][j])
                               ),
                        atol
                    ))
                {
                    return false;
                }
            }
        }
        return true;
    }

    // Getters
    std::vector<std::vector<double>> get_box_matrix() const
    {
        return _box_matrix;
    }
    std::vector<double> get_box_lengths() const { return _box_lengths; }
    std::vector<double> get_box_angles() const { return _box_angles; }
    double              get_x() const { return _box_lengths[0]; }
    double              get_y() const { return _box_lengths[1]; }
    double              get_z() const { return _box_lengths[2]; }
    double              get_alpha() const { return _box_angles[0]; }
    double              get_beta() const { return _box_angles[1]; }
    double              get_gamma() const { return _box_angles[2]; }

    // Setters
    void set_box_matrix(std::vector<std::vector<double>> box_matrix)
    {
        _box_matrix = box_matrix;
    }
    void set_box_lengths(std::vector<double> box_lengths)
    {
        _box_lengths = box_lengths;
    }
    void set_box_angles(std::vector<double> box_angles)
    {
        _box_angles = box_angles;
    }
    void set_x(double x) { _box_lengths[0] = x; }
    void set_y(double y) { _box_lengths[1] = y; }
    void set_z(double z) { _box_lengths[2] = z; }
    void set_alpha(double alpha) { _box_angles[0] = alpha; }
    void set_beta(double beta) { _box_angles[1] = beta; }
    void set_gamma(double gamma) { _box_angles[2] = gamma; }
};

#endif   // CELL_HPP