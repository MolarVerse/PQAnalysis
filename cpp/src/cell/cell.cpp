#include "cell.hpp"

using namespace std;
constexpr float M_PIf = (float) M_PI;

Cell::Cell(float a, float b, float c, float alpha, float beta, float gamma)
{
    _box_lengths = vector({a, b, c});
    _box_angles  = vector({alpha, beta, gamma});
    _box_matrix  = _setup_box_matrix();
}

vector<float> Cell::_setup_box_matrix()
{
    // Initialize box matrix with 3 rows and 3 columns
    vector<float> box_matrix(3 * 3, 0.0);

    // Calculate cosines and sines of angles
    float alpha = _box_angles[0], beta = _box_angles[1], gamma = _box_angles[2];
    float cos_alpha = cosf(alpha / 180.0f * M_PIf);
    float cos_beta  = cosf(beta / 180.0f * M_PIf);
    float cos_gamma = cosf(gamma / 180.0f * M_PIf);
    float sin_gamma = sinf(gamma / 180.0f * M_PIf);
    float sin_beta  = sinf(beta / 180.0f * M_PIf);

    float a = _box_lengths[0], b = _box_lengths[1], c = _box_lengths[2];

    // Assign values to box matrix
    box_matrix[(0 * 3) + 0] = a;
    box_matrix[(0 * 3) + 1] = b * cos_gamma;
    box_matrix[(0 * 3) + 2] = c * cos_beta;
    box_matrix[(1 * 3) + 1] = b * sin_gamma;
    box_matrix[(1 * 3) + 2] =
        c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    box_matrix[(2 * 3) + 2] =
        c * sqrt(
                sin_beta * sin_beta - (cos_alpha - cos_beta * cos_gamma) *
                                          (cos_alpha - cos_beta * cos_gamma) /
                                          (sin_gamma * sin_gamma)
            );

    return box_matrix;
}

vector<float> Cell::bounding_edges()
{
    // Initialize edges with 8 rows and 3 columns
    vector<float> edges(8 * 3, 0.0);

    // Values to iterate over
    float values[2] = {-0.5, 0.5};

    // Triple nested loop for x, y, z coordinates
    for (int i = 0; i < 2; i++)
    {
        float x = values[i];
        for (int j = 0; j < 2; j++)
        {
            float y = values[j];
            for (int k = 0; k < 2; k++)
            {
                float z = values[k];

                int idx = i * (4 * 3) + j * (2 * 3) + k;

                // Perform matrix multiplication manually
                for (int col = 0; col < 3; col++)
                {
                    edges[(idx * 3) + col] = _box_matrix[(0 * 3) + col] * x +
                                             _box_matrix[(1 * 3) + col] * y +
                                             _box_matrix[(2 * 3) + col] * z;
                }
            }
        }
    }

    return edges;
}

float Cell::volume() const
{
    // Calculate volume using determinant of box matrix
    return _box_matrix[(0 * 3) + 0] * _box_matrix[(1 * 3) + 1] *
               _box_matrix[(2 * 3) + 2] +
           _box_matrix[(0 * 3) + 1] * _box_matrix[(1 * 3) + 2] *
               _box_matrix[(2 * 3) + 0] +
           _box_matrix[(0 * 3) + 2] * _box_matrix[(1 * 3) + 0] *
               _box_matrix[(2 * 3) + 1] -
           _box_matrix[(0 * 3) + 2] * _box_matrix[(1 * 3) + 1] *
               _box_matrix[(2 * 3) + 0] -
           _box_matrix[(0 * 3) + 1] * _box_matrix[(1 * 3) + 0] *
               _box_matrix[(2 * 3) + 2] -
           _box_matrix[(0 * 3) + 0] * _box_matrix[(1 * 3) + 2] *
               _box_matrix[(2 * 3) + 1];
}

bool Cell::is_vacuum() const
{   // Check if volume is close to maximum float value
    return volume() >
           std::numeric_limits<float>::max() * 0.99 * 0.99 * 0.99 * 0.99;
}

vector<float> Cell::image(vector<float> pos)
{
    vector<float> image(pos.size(), 0.0);

    // Check for orthorhombic cell optimization
    bool is_orthorhombic =
        (_box_angles[0] == 90.0 && _box_angles[1] == 90.0 &&
         _box_angles[2] == 90.0);

    if (is_orthorhombic)
    {
        // Fast path for orthorhombic cells
        for (size_t i = 0; i < pos.size() / 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                float scaled = pos[(i * 3) + j] / _box_lengths[j];
                image[(i * 3) + j] =
                    pos[(i * 3) + j] - _box_lengths[j] * round(scaled);
            }
        }
    }
    else
    {
        // General case using box matrix
        vector<float> fractional(pos.size(), 0.0);

        // Convert to fractional coordinates
        for (size_t i = 0; i < pos.size() / 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    fractional[(i * 3) + j] +=
                        pos[(i * 3) + k] * _box_matrix[(k * 3) + j];
                }
            }
            // Wrap to [-0.5, 0.5)
            for (int j = 0; j < 3; j++)
            {
                fractional[(i * 3) + j] -= round(fractional[(i * 3) + j]);
            }
        }

        // Convert back to cartesian
        for (size_t i = 0; i < pos.size() / 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    image[(i * 3) + j] +=
                        fractional[(i * 3) + k] * _box_matrix[(j * 3) + k];
                }
            }
        }
    }

    return image;
}

Cell &Cell::init_from_box_matrix(vector<float> box_matrix)
{
    // Assign box matrix
    _box_matrix = box_matrix;

    // Calculate box lengths and angles
    _box_lengths = {
        sqrt(
            _box_matrix[(0 * 3) + 0] * _box_matrix[(0 * 3) + 0] +
            _box_matrix[(1 * 3) + 0] * _box_matrix[(1 * 3) + 0] +
            _box_matrix[(2 * 3) + 0] * _box_matrix[(2 * 3) + 0]
        ),
        sqrt(
            _box_matrix[(0 * 3) + 1] * _box_matrix[(0 * 3) + 1] +
            _box_matrix[(1 * 3) + 1] * _box_matrix[(1 * 3) + 1] +
            _box_matrix[(2 * 3) + 1] * _box_matrix[(2 * 3) + 1]
        ),
        sqrt(
            _box_matrix[(0 * 3) + 2] * _box_matrix[(0 * 3) + 2] +
            _box_matrix[(1 * 3) + 2] * _box_matrix[(1 * 3) + 2] +
            _box_matrix[(2 * 3) + 2] * _box_matrix[(2 * 3) + 2]
        )
    };

    _box_angles = {
        acos(
            (_box_matrix[(0 * 3) + 1] * _box_matrix[(0 * 3) + 2] +
             _box_matrix[(1 * 3) + 1] * _box_matrix[(1 * 3) + 2] +
             _box_matrix[(2 * 3) + 1] * _box_matrix[(2 * 3) + 2]) /
            (_box_lengths[1] * _box_lengths[2])
        ) * 180.0f /
            M_PIf,
        acos(
            (_box_matrix[(0 * 3) + 0] * _box_matrix[(0 * 3) + 2] +
             _box_matrix[(1 * 3) + 0] * _box_matrix[(1 * 3) + 2] +
             _box_matrix[(2 * 3) + 0] * _box_matrix[(2 * 3) + 2]) /
            (_box_lengths[0] * _box_lengths[2])
        ) * 180.0f /
            M_PIf,
        acos(
            (_box_matrix[(0 * 3) + 0] * _box_matrix[(0 * 3) + 1] +
             _box_matrix[(1 * 3) + 0] * _box_matrix[(1 * 3) + 1] +
             _box_matrix[(2 * 3) + 0] * _box_matrix[(2 * 3) + 1]) /
            (_box_lengths[0] * _box_lengths[1])
        ) * 180.0f /
            M_PIf
    };

    return *this;
}