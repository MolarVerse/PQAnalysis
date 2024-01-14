"""
A module containing different decorators associated to methods of the AtomicSystem class.

...

Functions
---------
check_atoms_pos
    Decorator which checks that the number of atoms is equal to the number of positions.
check_atoms_has_mass
    Decorator which checks that all atoms have mass information.
"""

from decorator import decorator

from ..exceptions import AtomicSystemPositionsError, AtomicSystemMassError


@decorator
def check_atoms_pos(func, *args, **kwargs):
    """
    Decorator which checks that the number of atoms is equal to the number of positions.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    AtomicSystemPositionsError
        If the number of atoms is not equal the number of positions.
    """
    self = args[0]

    if self.pos.shape[0] != len(self.atoms):
        raise AtomicSystemPositionsError()

    return func(*args, **kwargs)


@decorator
def check_atoms_has_mass(func, *args, **kwargs):
    """
    Decorator which checks that all atoms have mass information.

    Parameters
    ----------
    func : function
        Function to be decorated.

    Raises
    ------
    ValueError
        If any atom does not have mass information.
    """

    self = args[0]

    if not all([atom.mass is not None for atom in self.atoms]):
        raise AtomicSystemMassError()

    return func(*args, **kwargs)
