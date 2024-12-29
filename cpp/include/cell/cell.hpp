#ifndef CELL_HPP
#define CELL_HPP

#include <math.h>

#include <string>
#include <vector>

class Cell
{
   private:
    std::vector<float> _box_lengths, _box_angles;
    std::vector<float> _box_matrix;

   public:
    Cell(float x, float y, float z, float alpha, float beta, float gamma);

   private:
    std::vector<float> _setup_box_matrix();

   public:
    std::vector<float> bounding_edges();
    float              volume() const;
    bool               is_vacuum() const;
    std::vector<float> image(std::vector<float> pos);
    Cell              &init_from_box_matrix(std::vector<float> box_matrix);

    bool isclose(const Cell &other, float rtol = 1e-9, float atol = 0.0) const
    {
        if (is_vacuum() && other.is_vacuum())
        {
            return true;
        }

        // Compare box matrix with other box matrix contents
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (fabs(_box_matrix[i + j] - other._box_matrix[i + j]) >
                    fmax(
                        rtol * fmax(
                                   fabs(_box_matrix[i + j]),
                                   fabs(other._box_matrix[i + j])
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
    std::vector<float> get_box_matrix() const { return _box_matrix; }
    std::vector<float> get_box_lengths() const { return _box_lengths; }
    std::vector<float> get_box_angles() const { return _box_angles; }
    float              get_x() const { return _box_lengths[0]; }
    float              get_y() const { return _box_lengths[1]; }
    float              get_z() const { return _box_lengths[2]; }
    float              get_alpha() const { return _box_angles[0]; }
    float              get_beta() const { return _box_angles[1]; }
    float              get_gamma() const { return _box_angles[2]; }

    // Setters
    void set_box_matrix(std::vector<float> box_matrix)
    {
        _box_matrix = box_matrix;
    }
    void set_box_lengths(std::vector<float> box_lengths)
    {
        _box_lengths = box_lengths;
    }
    void set_box_angles(std::vector<float> box_angles)
    {
        _box_angles = box_angles;
    }
    void set_x(float x) { _box_lengths[0] = x; }
    void set_y(float y) { _box_lengths[1] = y; }
    void set_z(float z) { _box_lengths[2] = z; }
    void set_alpha(float alpha) { _box_angles[0] = alpha; }
    void set_beta(float beta) { _box_angles[1] = beta; }
    void set_gamma(float gamma) { _box_angles[2] = gamma; }

    // string representation
    std::string to_string() const
    {
        if (is_vacuum())
        {
            return "Cell()";
        }

        return "Cell(" + std::to_string(_box_lengths[0]) + ", " +
               std::to_string(_box_lengths[1]) + ", " +
               std::to_string(_box_lengths[2]) + ", " +
               std::to_string(_box_angles[0]) + ", " +
               std::to_string(_box_angles[1]) + ", " +
               std::to_string(_box_angles[2]) + ")";
    }
};

#endif   // CELL_HPP