import numpy as np

from PQAnalysis.utils.math import allclose_vectorized

from . import pytestmark



def test_allclose_vectorized():
    a = np.array([1, 2, 3])
    b = np.array([1, 2, 3])

    assert allclose_vectorized(a, b)

    a = np.array([1, 2, 3])
    b = np.array([1, 2, 3.0000009])

    assert not allclose_vectorized(a, b)
    assert allclose_vectorized(a, b, atol=1e-6)
    assert not allclose_vectorized(a, b, atol=1e-7)
    assert allclose_vectorized(a, b, rtol=1e-6)
    assert not allclose_vectorized(a, b, rtol=1e-7)

    a = np.array([1, 2, 10])
    b = np.array([1, 2, 10.000001])

    assert not allclose_vectorized(a, b, atol=0.9e-6)
    assert allclose_vectorized(a, b, atol=1e-6)
    assert allclose_vectorized(a, b, rtol=1e-6)
    assert allclose_vectorized(a, b, rtol=1.1e-6)
    assert not allclose_vectorized(a, b, rtol=0.9e-7)

    assert allclose_vectorized(a, b, atol=0.9e-6, rtol=1e-6)
    assert allclose_vectorized(a, b, atol=1e-6, rtol=0.9e-7)
    assert not allclose_vectorized(a, b, atol=0.9e-6, rtol=0.9e-7)



# Compare this snippet from tests/utils/test_math.py:
