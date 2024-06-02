"""
A module containing type hints for the PQAnalysis package.
"""
from numbers import Real
from typing import Annotated

import numpy as np

from beartype.vale import Is
from beartype.typing import NewType

#: A type hint for positive integers
PositiveInt = NewType(
    "PositiveInt",
    Annotated[int | np.int_, Is[lambda int: int > 0]],
)

# :A type hint for positive real numbers
PositiveReal = NewType(
    "PositiveReal", Annotated[
        Real,
        Is[lambda real: real >= 0.0],
    ]
)

#: A type variable for a 1D np.ndarray with dtype np.number
Np1DNumberArray = NewType(
    "Np1DNumberArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim == 1 and
           (np.issubdtype(array.dtype, np.number)) or len(array) == 0],
    ]
)

#: A type hint for a 2D np.ndarray with dtype np.number
Np2DNumberArray = NewType(
    "Np2DNumberArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim == 2 and
           (np.issubdtype(array.dtype, np.number))],
    ]
)

#: A type hint for a nD np.ndarray with dtype np.number
NpnDNumberArray = NewType(
    "NpnDNumberArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim > 0 and
           (np.issubdtype(array.dtype, np.number))],
    ]
)

#: A type hint for a 3x3 np.ndarray matrix with dtype np.number
Np3x3NumberArray = NewType(
    "Np2x2NumberArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim == 2 and array.shape == (3, 3) and
           (np.issubdtype(array.dtype, np.number))],
    ]
)

#: A type hint for a 2D np.ndarray with dtype np.integer
Np2DIntArray = NewType(
    "Np2DIntArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim == 2 and
           (np.issubdtype(array.dtype, np.integer))],
    ]
)

#: A type hint for a 1D np.ndarray with dtype np.integer
Np1DIntArray = NewType(
    "Np1DIntArray",
    Annotated[
        np.ndarray,
        Is[lambda array: array.ndim == 1 and
           (np.issubdtype(array.dtype, np.integer)) or len(array) == 0],
    ]
)

#: A type hint for a range object
Range = NewType("Range", Annotated[range, Is[lambda r: isinstance(r, range)]])

Bool = NewType(
    "Bool",
    Annotated[
        bool | np.bool_,
        Is[lambda b: isinstance(b, (bool, np.bool_))],
    ]
)
