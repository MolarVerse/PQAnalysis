"""
A module containing different unit conversions and constants.
"""

from numbers import Real

from PQAnalysis.exceptions import PQValueError

mol = 6.02214076e23  # pylint: disable=invalid-name

J = 1.0
J_per_mol = J * mol  # pylint: disable=invalid-name
cal = J / 4.184  # pylint: disable=invalid-name
kcal = cal / 1000.0  # pylint: disable=invalid-name
eV = 6.241506363094e18 * J  # pylint: disable=invalid-name
kcal_per_mol = kcal * mol  # pylint: disable=invalid-name


def calculate_simulation_time(
    timesteps: Real,
    timestep_fs: Real,
    unit: str = "ps",
) -> float:
    """
    Calculates simulation time from a timestep count.

    Parameters
    ----------
    timesteps : Real
        The number of timesteps.
    timestep_fs : Real
        The timestep size in femtoseconds.
    unit : str, optional
        The output unit. Supported values are ``"fs"`` and ``"ps"``, by
        default ``"ps"``.

    Returns
    -------
    float
        The simulation time in the requested unit.

    Raises
    ------
    PQValueError
        If timesteps or timestep_fs are negative, or if unit is unsupported.
    """
    if timesteps < 0:
        raise PQValueError("timesteps has to be greater than or equal to 0.")

    if timestep_fs < 0:
        raise PQValueError("timestep_fs has to be greater than or equal to 0.")

    time_fs = float(timesteps * timestep_fs)

    if unit.lower() == "fs":
        return time_fs

    if unit.lower() == "ps":
        return time_fs / 1000.0

    raise PQValueError(
        f"Unsupported time unit {unit}. Supported units are fs and ps."
    )
