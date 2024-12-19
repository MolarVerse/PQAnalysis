#include "cell.hpp"

using namespace std;

CoreCell::CoreCell()
{
    _box_lengths = array_d::zeros(3);
    _box_angles = array_d::zeros(3);
    _box_matrix = Eigen::Matrix::Zero(3, 3);
}

CoreCell::CoreCell(double a, double b, double c, double alpha, double beta, double gamma)
{
    _box_lengths = {a, b, c};
    _box_angles = {alpha, beta, gamma};
    _box_matrix = _setup_box_matrix();
}

Eigen::Matrix CoreCell::_setup_box_matrix()
{
    // Initialize box matrix with 3 rows and 3 columns
    Eigen::Matrix box_matrix = Eigen::Matrix::Zero(3, 3);

    // Calculate cosines and sines of angles
    double alpha = _box_angles[0], beta = _box_angles[1], gamma = _box_angles[2];
    double cos_alpha = cos(alpha / 180.0 * M_PI);
    double cos_beta = cos(beta / 180.0 * M_PI);
    double cos_gamma = cos(gamma / 180.0 * M_PI);
    double sin_gamma = sin(gamma / 180.0 * M_PI);
    double sin_beta = sin(beta / 180.0 * M_PI);

    double a = _box_lengths[0], b = _box_lengths[1], c = _box_lengths[2];

    // Assign values to box matrix
    box_matrix[0][0] = a;
    box_matrix[0][1] = b * cos_gamma;
    box_matrix[0][2] = c * cos_beta;
    box_matrix[1][1] = b * sin_gamma;
    box_matrix[1][2] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    box_matrix[2][2] = c * sqrt(
                               sin_beta * sin_beta -
                               (cos_alpha - cos_beta * cos_gamma) * (cos_alpha - cos_beta * cos_gamma) / (sin_gamma * sin_gamma));

    return box_matrix;
}

Eigen::Matrix CoreCell::bouding_edges()
{
    // Initialize edges array with 8 rows and 3 columns - array_t<double> is a 2D array
    Matrix edges = Matrix::Zero(8, 3);

    // Values to iterate over
    vector_d values = {-0.5, 0.5};
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
                array_d vec = {x, y, z};
                Matrix result = _box_matrix * vec;

                // Assign result to edges
                edges.row(idx) = result;
            }
        }
    }

    return edges;
}

double CoreCell::volume()
{
    // Calculate volume using determinant of box matrix
    return _box_matrix.determinant();
}

bool CoreCell::is_vacuum()
{
    // Check if volume is zero
    return volume() == 0;
}

array_d CoreCell::image(array_d pos)
{

    // Unpack position array
    auto position = array_d::ensure(pos);

    // Initialize image matrix with same size as position
    Matrix image = Eigen::Matrix::Zero(position.shape(0), position.shape(1));

    // Iterate over all positions
    for (int i = 0; i < position.shape(0); i++)
    {
        // Perform matrix multiplication manually
        for (int col = 0; col < position.shape(1); col++)
        {
            image[i][col] = _box_matrix[col][0] * position(i, 0) +
                            _box_matrix[col][1] * position(i, 1) +
                            _box_matrix[col][2] * position(i, 2);
        }
    }

    return image;
}

CoreCell &CoreCell::init_from_box_matrix(array_d box_matrix)
{
    // Convert numpy array to Eigen matrix
    _box_matrix = numpy_to_eigen(box_matrix);

    // Calculate box lengths and angles
    _box_lengths = {sqrt(_box_matrix.col(0).squaredNorm()),
                    sqrt(_box_matrix.col(1).squaredNorm()),
                    sqrt(_box_matrix.col(2).squaredNorm())};

    _box_angles = {acos(_box_matrix.col(1).dot(_box_matrix.col(2)) / (_box_lengths[1] * _box_lengths[2])) * 180.0 / M_PI,
                   acos(_box_matrix.col(0).dot(_box_matrix.col(2)) / (_box_lengths[0] * _box_lengths[2])) * 180.0 / M_PI,
                   acos(_box_matrix.col(0).dot(_box_matrix.col(1)) / (_box_lengths[0] * _box_lengths[1])) * 180.0 / M_PI};

    return this;
}