#include "cell.hpp"

using namespace std;

Cell::Cell()
{
    _box_lengths = vector(3, 0.0);
    _box_angles  = vector(3, 0.0);
    _box_matrix  = vector(3, vector(3, 0.0));
}

Cell::Cell(double x, double y, double z)
{
    _box_lengths = vector({x, y, z});
    _box_angles  = vector(3, 90.0);
    _box_matrix  = _setup_box_matrix();
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

array_d Cell::bouding_edges()
{
    // Initialize edges array with 8 rows and 3 columns - array_t<double> is a
    // 2D array
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

                // Calculate index similar to Python version
                int idx = i * 4 + j * 2 + k;

                // Perform matrix multiplication using Eigen
                // Perform matrix multiplication using Eigen
                vector<double> result(3);
                result[0] = _box_matrix[0][0] * x + _box_matrix[0][1] * y +
                            _box_matrix[0][2] * z;
                result[1] = _box_matrix[1][0] * x + _box_matrix[1][1] * y +
                            _box_matrix[1][2] * z;
                result[2] = _box_matrix[2][0] * x + _box_matrix[2][1] * y +
                            _box_matrix[2][2] * z;

                // Assign result to edges
                edges[idx] = result;
            }
        }
    }

    return vector_2d_to_array(edges);
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

array_d Cell::image(array_d pos)
{
    // Unpack position array
    auto position = array_to_vector_2d(pos);

    // Initialize image matrix with same size as position
    vector<vector<double>> image(position.size(), vector<double>(3));

    // Iterate over all positions
    for (int i = 0; i < position.size(); i++)
    {
        // Perform matrix multiplication manually
        for (int col = 0; col < 3; col++)
        {
            image[i][col] = _box_matrix[0][col] * position[i][0] +
                            _box_matrix[1][col] * position[i][1] +
                            _box_matrix[2][col] * position[i][2];
        }
    }

    return vector_2d_to_array(image);
}

Cell &Cell::init_from_box_matrix(array_d box_matrix)
{
    // Assign box matrix
    _box_matrix = array_to_vector_2d(box_matrix);

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

array_d vector_2d_to_array(const std::vector<std::vector<double>> &vec)
{
    // Determine the dimensions of the 2D array
    size_t rows = vec.size();
    size_t cols = vec.empty() ? 0 : vec[0].size();

    // Flatten the 2D vector into a 1D vector
    std::vector<double> flattened;
    flattened.reserve(rows * cols);
    for (const auto &row : vec)
    {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }

    // Create a NumPy array from the flattened data
    array_d array({rows, cols}, flattened.data());
    return array;
}

std::vector<std::vector<double>> array_to_vector_2d(array_d arr)
{
    // Get the dimensions of the NumPy array
    py::buffer_info info = arr.request();
    size_t          rows = info.shape[0];
    size_t          cols = info.shape[1];

    // Create a 2D vector from the NumPy array
    std::vector<std::vector<double>> vec(rows, std::vector<double>(cols));
    double                          *ptr = static_cast<double *>(info.ptr);
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            vec[i][j] = ptr[i * cols + j];
        }
    }

    return vec;
}

std::vector<double> array_to_vector_1d(array_d arr)
{
    // Get the dimensions of the NumPy array
    py::buffer_info info = arr.request();
    size_t          size = info.size;

    // Create a 1D vector from the NumPy array
    std::vector<double> vec(size);
    double             *ptr = static_cast<double *>(info.ptr);
    for (size_t i = 0; i < size; i++)
    {
        vec[i] = ptr[i];
    }

    return vec;
}