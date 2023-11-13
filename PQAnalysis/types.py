import numpy as np

from beartype.vale import Is
from typing import Annotated


Np2DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 2 and
                                           (np.issubdtype(array.dtype, np.number))]]

Np1DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 1 and
                                           (np.issubdtype(array.dtype, np.number))]]

Np3x3NumberArray = Annotated[np.ndarray, Is[lambda array:
                                            array.ndim == 2 and
                                            array.shape == (3, 3) and
                                            (np.issubdtype(array.dtype, np.number))]]

Np2DIntArray = Annotated[np.ndarray, Is[lambda array:
                                        array.ndim == 2 and
                                        (np.issubdtype(array.dtype, np.integer))]]
