"""
A collection of utility functions for mathematical operations.
"""

import math
import numpy as np

from PQAnalysis.types import PositiveReal, Bool



def allclose_vectorized(
    a,
    b,
    rtol: PositiveReal = 1e-09,
    atol: PositiveReal = 0.0,
) -> Bool:
    """
    Perform element-wise comparison of two arrays to determine if they
    are equal within a tolerance. This function is a vectorized version
    of `math.isclose` that can be applied to numpy arrays.
    
    This function uses the full symmetric definition of `math.isclose` to
    compare two arrays, by applying the following formula element-wise:
    
    .. math::
        
            |a - b| <= max(rtol * max(|a|, |b|), atol))
    
    For a full discussion on why this function is needed, and replaces `np.allclose`,
    see [pull request #74](https://github.com/MolarVerse/PQAnalysis/pull/74)

    Parameters
    ----------
    a: np.ndarray
        first array to compare
    b : np.ndarray
        second array to compare
    rtol : PositiveReal, optional
        the relative tolerance parameter, by default 1e-09
    atol : PositiveReal, optional
        the absolute tolerance parameter, by default 0.0

    Returns
    -------
    Bool
        True if the arrays are element-wise equal within the given tolerance, False otherwise
    """

    # otypes needed: ValueError: cannot call `vectorize` on size 0 inputs unless `otypes` is set
    isclose_vectorized = np.vectorize(math.isclose, otypes=[bool])

    result = isclose_vectorized(a, b, rel_tol=rtol, abs_tol=atol)

    return np.all(result)
