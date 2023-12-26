import numpy as np
import _io

from numbers import Real
from beartype.vale import Is
from typing import Annotated

PositiveInt = Annotated[int, Is[lambda int: int > 0]]

PositiveReal = Annotated[Real, Is[lambda real: real >= 0.0]]

Np2DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 2 and
                                           (np.issubdtype(array.dtype, np.number))]]

Np1DNumberArray = Annotated[np.ndarray, Is[lambda array:
                                           array.ndim == 1 and
                                           (np.issubdtype(array.dtype, np.number)) or
                                           len(array) == 0]]

Np3x3NumberArray = Annotated[np.ndarray, Is[lambda array:
                                            array.ndim == 2 and
                                            array.shape == (3, 3) and
                                            (np.issubdtype(array.dtype, np.number))]]

Np2DIntArray = Annotated[np.ndarray, Is[lambda array:
                                        array.ndim == 2 and
                                        (np.issubdtype(array.dtype, np.integer))]]

Np1DIntArray = Annotated[np.ndarray, Is[lambda array:
                                        array.ndim == 1 and
                                        (np.issubdtype(array.dtype, np.integer)) or
                                        len(array) == 0]]

FILE = _io.TextIOWrapper

Range = Annotated[range, Is[lambda r: isinstance(r, range)]]
Bool = Annotated[bool, Is[lambda b: isinstance(b, bool)]]
PositiveInt = Annotated[int, Is[lambda i: isinstance(i, int) and i > 0]]
