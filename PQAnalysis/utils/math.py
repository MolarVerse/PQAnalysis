import numpy as np


def to_numpy_array(x):
    return np.asarray([x]) if np.isscalar(x) else np.asarray(x)
