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