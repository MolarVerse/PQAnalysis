#include "cell.hpp"

using namespace std;

Cell::Cell()
{
    _box_lengths = vector(3, 0.0);
    _box_angles  = vector(3, 0.0);
    _box_matrix  = vector(3, vector(3, 0.0));
}

Cell::Cell(
    double a,
    double b,
    double c,
    double alpha,
    double beta,
    double gamma
)
{
    _box_lengths = vector({a, b, c});
    _box_angles  = vector({alpha, beta, gamma});
    _box_matrix  = _setup_box_matrix();
}

vector<vector<double>> Cell::_setup_box_matrix()
{
    // Initialize box matrix with 3 rows and 3 columns
    vector<vector<double>> box_matrix(3, vector<double>(3));

    // Calculate cosines and sines of angles
    double alpha = _box_angles[0], beta = _box_angles[1],
           gamma     = _box_angles[2];
    double cos_alpha = cos(alpha / 180.0 * M_PI);
    double cos_beta  = cos(beta / 180.0 * M_PI);
    double cos_gamma = cos(gamma / 180.0 * M_PI);
    double sin_gamma = sin(gamma / 180.0 * M_PI);
    double sin_beta  = sin(beta / 180.0 * M_PI);

    double a = _box_lengths[0], b = _box_lengths[1], c = _box_lengths[2];

    // Assign values to box matrix
    box_matrix[0][0] = a;
    box_matrix[0][1] = b * cos_gamma;
    box_matrix[0][2] = c * cos_beta;
    box_matrix[1][1] = b * sin_gamma;
    box_matrix[1][2] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    box_matrix[2][2] =
        c * sqrt(
                sin_beta * sin_beta - (cos_alpha - cos_beta * cos_gamma) *
                                          (cos_alpha - cos_beta * cos_gamma) /
                                          (sin_gamma * sin_gamma)
            );

    return box_matrix;
}

vector<vector<double>> Cell::bounding_edges()
{
    // Initialize edges with 8 rows and 3 columns
    vector<vector<double>> edges(8, vector<double>(3));

    // Values to iterate over
    double values[2] = {-0.5, 0.5};

    // Triple nested loop for x, y, z coordinates
    for (int i = 0; i < 2; i++)
    {
        double x = values[i];
        for (int j = 0; j < 2; j++)
        {
            double y = values[j];
            for (int k = 0; k < 2; k++)
            {
                double z = values[k];

                int idx = i * 4 + j * 2 + k;

                // Perform matrix multiplication manually
                for (int col = 0; col < 3; col++)
                {
                    edges[idx][col] = _box_matrix[0][col] * x +
                                      _box_matrix[1][col] * y +
                                      _box_matrix[2][col] * z;
                }
            }
        }
    }

    return edges;
}

double Cell::volume()
{
    // Calculate volume using determinant of box matrix
    return _box_matrix[0][0] * _box_matrix[1][1] * _box_matrix[2][2] +
           _box_matrix[0][1] * _box_matrix[1][2] * _box_matrix[2][0] +
           _box_matrix[0][2] * _box_matrix[1][0] * _box_matrix[2][1] -
           _box_matrix[0][2] * _box_matrix[1][1] * _box_matrix[2][0] -
           _box_matrix[0][1] * _box_matrix[1][0] * _box_matrix[2][2] -
           _box_matrix[0][0] * _box_matrix[1][2] * _box_matrix[2][1];
}

bool Cell::is_vacuum()
{
    // Check if volume is zero
    return volume() == 0;
}

vector<vector<double>> Cell::image(vector<vector<double>> pos)
{
    vector<vector<double>> image(pos.size(), vector<double>(3, 0.0));

    // Check for orthorhombic cell optimization
    bool is_orthorhombic =
        (_box_angles[0] == 90.0 && _box_angles[1] == 90.0 &&
         _box_angles[2] == 90.0);

    if (is_orthorhombic)
    {
        // Fast path for orthorhombic cells
        for (size_t i = 0; i < pos.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double scaled = pos[i][j] / _box_lengths[j];
                image[i][j]   = pos[i][j] - _box_lengths[j] * round(scaled);
            }
        }
    }
    else
    {
        // General case using box matrix
        vector<vector<double>> fractional(pos.size(), vector<double>(3, 0.0));

        // Convert to fractional coordinates
        for (size_t i = 0; i < pos.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    fractional[i][j] += pos[i][k] * _box_matrix[k][j];
                }
            }
            // Wrap to [-0.5, 0.5)
            for (int j = 0; j < 3; j++)
            {
                fractional[i][j] -= round(fractional[i][j]);
            }
        }

        // Convert back to cartesian
        for (size_t i = 0; i < pos.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    image[i][j] += fractional[i][k] * _box_matrix[j][k];
                }
            }
        }
    }

    return image;
}

Cell &Cell::init_from_box_matrix(vector<vector<double>> box_matrix)
{
    // Assign box matrix
    _box_matrix = box_matrix;

    // Calculate box lengths and angles
    _box_lengths = {
        sqrt(
            _box_matrix[0][0] * _box_matrix[0][0] +
            _box_matrix[1][0] * _box_matrix[1][0] +
            _box_matrix[2][0] * _box_matrix[2][0]
        ),
        sqrt(
            _box_matrix[0][1] * _box_matrix[0][1] +
            _box_matrix[1][1] * _box_matrix[1][1] +
            _box_matrix[2][1] * _box_matrix[2][1]
        ),
        sqrt(
            _box_matrix[0][2] * _box_matrix[0][2] +
            _box_matrix[1][2] * _box_matrix[1][2] +
            _box_matrix[2][2] * _box_matrix[2][2]
        )
    };

    _box_angles = {
        acos(
            (_box_matrix[0][1] * _box_matrix[0][2] +
             _box_matrix[1][1] * _box_matrix[1][2] +
             _box_matrix[2][1] * _box_matrix[2][2]) /
            (_box_lengths[1] * _box_lengths[2])
        ) * 180.0 /
            M_PI,
        acos(
            (_box_matrix[0][0] * _box_matrix[0][2] +
             _box_matrix[1][0] * _box_matrix[1][2] +
             _box_matrix[2][0] * _box_matrix[2][2]) /
            (_box_lengths[0] * _box_lengths[2])
        ) * 180.0 /
            M_PI,
        acos(
            (_box_matrix[0][0] * _box_matrix[0][1] +
             _box_matrix[1][0] * _box_matrix[1][1] +
             _box_matrix[2][0] * _box_matrix[2][1]) /
            (_box_lengths[0] * _box_lengths[1])
        ) * 180.0 /
            M_PI
    };

    return *this;
}