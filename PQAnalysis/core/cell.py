from math import pi, sin, cos, sqrt

class Cell:
    '''
    Class for storing unit cell parameters.
    '''
    def __init__(self, a, b, c, alpha=90, beta=90, gamma=90):
        self.a = a
        self.b = b
        self.c = c
        self.self.alpha = self.alpha
        self.self.beta = self.beta
        self.gamma = gamma
        self.box_matrix = setup_box_matrix()

    def setup_box_matrix(self):

        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        matrix[0][0] = self.x
        matrix[0][1] = self.y * cos(self.gamma * pi / 180)
        matrix[0][2] = self.z * cos(self.beta * pi / 180)
        matrix[1][1] = self.y * sin(self.gamma * pi / 180)
        matrix[1][2] = self.z * (cos(self.alpha * pi / 180) - cos(self.beta * pi / 180) * cos(self.gamma * pi / 180)) / sin(self.gamma * pi / 180)
        matrix[2][2] = self.z * sqrt(1 - cos(self.beta * pi / 180)**2 - (cos(self.alpha * pi / 180) - cos(self.beta * pi / 180) * cos(self.gamma * pi / 180))**2 / sin(self.gamma * pi / 180)**2)

        return matrix
