import numpy as np


class Cell:
    '''
    Class for storing unit cell parameters.
    '''

    def __init__(self, x, y, z, alpha=90, beta=90, gamma=90):
        self.x = x
        self.y = y
        self.z = z
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.box_matrix = self.setup_box_matrix()

    def setup_box_matrix(self):

        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        alpha = np.deg2rad(self.alpha)
        beta = np.deg2rad(self.beta)
        gamma = np.deg2rad(self.gamma)

        matrix[0][0] = self.x
        matrix[0][1] = self.y * np.cos(gamma)
        matrix[0][2] = self.z * np.cos(beta)
        matrix[1][1] = self.y * np.sin(gamma)
        matrix[1][2] = self.z * \
            (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        matrix[2][2] = self.z * np.sqrt(1 - np.cos(beta)**2 - (
            np.cos(alpha) - np.cos(beta) * np.cos(gamma))**2 / np.sin(gamma)**2)

        return matrix

    @property
    def bounding_edges(self):
        edges = np.zeros((8, 3))
        for i, x in enumerate([-0.5, 0.5]):
            for j, y in enumerate([-0.5, 0.5]):
                for k, z in enumerate([-0.5, 0.5]):
                    edges[i*4+j*2+k, :] = self.box_matrix @ np.array([x, y, z])

        return edges

    @property
    def volume(self):
        return np.linalg.det(self.box_matrix)

    @property
    def box_lengths(self):
        return np.array([self.x, self.y, self.z])

    @property
    def box_angles(self):
        return np.array([self.alpha, self.beta, self.gamma])

    def image(self, pos):
        '''
        Returns the image of a position vector in the unit cell.
        '''
        if np.shape(pos) == (3,):
            pos = np.reshape(pos, (1, 3))

        fractional_pos = [np.linalg.inv(self.box_matrix) @ i for i in pos]

        fractional_pos -= np.round(fractional_pos)

        pos = [self.box_matrix @ i for i in fractional_pos]

        if np.shape(pos) == (1, 3):
            pos = np.reshape(pos, (3,))

        return pos
