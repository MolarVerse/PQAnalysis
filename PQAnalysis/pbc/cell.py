from math import pi, sin, cos, sqrt
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

        matrix[0][0] = self.x
        matrix[0][1] = self.y * cos(self.gamma * pi / 180)
        matrix[0][2] = self.z * cos(self.beta * pi / 180)
        matrix[1][1] = self.y * sin(self.gamma * pi / 180)
        matrix[1][2] = self.z * (cos(self.alpha * pi / 180) - cos(
            self.beta * pi / 180) * cos(self.gamma * pi / 180)) / sin(self.gamma * pi / 180)
        matrix[2][2] = self.z * sqrt(1 - cos(self.beta * pi / 180)**2 - (cos(self.alpha * pi / 180) - cos(
            self.beta * pi / 180) * cos(self.gamma * pi / 180))**2 / sin(self.gamma * pi / 180)**2)

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
