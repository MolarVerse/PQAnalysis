from math import pi, sin, cos, sqrt
import numpy as np


class Cell:
    '''
    Class for storing unit cell parameters.
    '''

    def __init__(self, a, b, c, alpha=90, beta=90, gamma=90):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.box_matrix = self.setup_box_matrix()

    def setup_box_matrix(self):

        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        matrix[0][0] = self.a
        matrix[0][1] = self.b * cos(self.gamma * pi / 180)
        matrix[0][2] = self.c * cos(self.beta * pi / 180)
        matrix[1][1] = self.b * sin(self.gamma * pi / 180)
        matrix[1][2] = self.c * (cos(self.alpha * pi / 180) - cos(
            self.beta * pi / 180) * cos(self.gamma * pi / 180)) / sin(self.gamma * pi / 180)
        matrix[2][2] = self.c * sqrt(1 - cos(self.beta * pi / 180)**2 - (cos(self.alpha * pi / 180) - cos(
            self.beta * pi / 180) * cos(self.gamma * pi / 180))**2 / sin(self.gamma * pi / 180)**2)

        return matrix
