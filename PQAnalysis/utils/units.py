"""
A module containing different units related to the unum subpackage.
"""

mol = 6.02214076e23  # pylint: disable=invalid-name

J = 1.0
J_per_mol = J * mol  # pylint: disable=invalid-name
cal = J / 4.184  # pylint: disable=invalid-name
kcal = cal / 1000.0  # pylint: disable=invalid-name
eV = 6.241506363094e18 * J  # pylint: disable=invalid-name
kcal_per_mol = kcal * mol  # pylint: disable=invalid-name
