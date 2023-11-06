import numpy as np

from beartype.vale import Is
from typing import Annotated

# Import the requisite machinery.
from beartype import beartype, BeartypeConf


Numpy2DFloatArray = Annotated[np.ndarray, Is[lambda array:
                                             array.ndim == 2 and
                                             (np.issubdtype(array.dtype, np.floating) or np.issubdtype(array.dtype, np.integer))]]

Numpy1DFloatArray = Annotated[np.ndarray, Is[lambda array:
                                             array.ndim == 1 and
                                             (np.issubdtype(array.dtype, np.floating) or np.issubdtype(array.dtype, np.integer))]]

Numpy3x3FloatArray = Annotated[np.ndarray, Is[lambda array:
                                              array.ndim == 2 and
                                              array.shape == (3, 3) and
                                              (np.issubdtype(array.dtype, np.floating) or np.issubdtype(array.dtype, np.integer))]]
