import numpy as np

from numbers import Real
from beartype.vale import Is
from typing import Annotated

"""
A type hint for positive integers
"""
PositiveInt = Annotated[int, Is[lambda int: int > 0]]

"""
A type hint for positive real numbers
"""
PositiveReal = Annotated[Real, Is[lambda real: real >= 0.0]]

"""
A type hint for a 2D np.ndarray with dtype np.number
"""
Np2DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 2 and
                                           (np.issubdtype(array.dtype, np.number))]]

"""
A type hint for a 1D np.ndarray with dtype np.number
"""
Np1DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 1 and
                                           (np.issubdtype(array.dtype, np.number)) or
                                           len(array) == 0]]

"""
A type hint for a 3x3 np.ndarray matrix with dtype np.number
"""
Np3x3NumberArray = Annotated[np.ndarray, Is[lambda array:
                                            array.ndim == 2 and
                                            array.shape == (3, 3) and
                                            (np.issubdtype(array.dtype, np.number))]]

"""
A type hint for a 2D np.ndarray with dtype np.integer
"""
Np2DIntArray = Annotated[np.ndarray, Is[lambda array:
                                        array.ndim == 2 and
                                        (np.issubdtype(array.dtype, np.integer))]]

"""
A type hint for a 1D np.ndarray with dtype np.integer
"""
Np1DIntArray = Annotated[np.ndarray, Is[lambda array:
                                        array.ndim == 1 and
                                        (np.issubdtype(array.dtype, np.integer)) or
                                        len(array) == 0]]

"""
A type hint for a range object
"""
Range = Annotated[range, Is[lambda r: isinstance(r, range)]]
