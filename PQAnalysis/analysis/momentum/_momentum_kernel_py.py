"""Python fallback for the legacy-compatible momentum kernel."""

import math


def _legacy_norm3(x, y, z):
    """Replicates Julia's ``generic_norm2`` for three float64 values."""
    max_abs = max(abs(x), abs(y), abs(z))

    if max_abs == 0.0 or math.isinf(max_abs) or math.isnan(max_abs):
        return max_abs

    square = max_abs * max_abs

    if math.isfinite((3.0 * max_abs) * max_abs) and square != 0.0:
        total = x * x
        total = total + y * y
        total = total + z * z
        return math.sqrt(total)

    value = abs(x) / max_abs
    total = value * value
    value = abs(y) / max_abs
    total = total + value * value
    value = abs(z) / max_abs
    total = total + value * value

    return max_abs * math.sqrt(total)


def legacy_momentum_norm(values, indices, masses, scale):
    """Calculates one scaled momentum norm in legacy operation order."""
    if len(masses) != len(indices):
        raise ValueError("indices and masses must have the same length")

    momentum_x = 0.0
    momentum_y = 0.0
    momentum_z = 0.0

    for index, mass_value in zip(indices, masses):
        row = int(index)
        mass = float(mass_value)
        momentum_x = momentum_x + mass * float(values[row, 0])
        momentum_y = momentum_y + mass * float(values[row, 1])
        momentum_z = momentum_z + mass * float(values[row, 2])

    return _legacy_norm3(
        momentum_x,
        momentum_y,
        momentum_z,
    ) * float(scale)
