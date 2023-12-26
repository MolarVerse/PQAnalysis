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

from .. import AtomicSystemPositionsError, AtomicSystemMassError


def check_atoms_pos(func):
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
    def wrapper(*args, **kwargs):

        if args[0].pos.shape[0] != len(args[0].atoms):
            raise AtomicSystemPositionsError()

        return func(*args, **kwargs)

    return wrapper


def check_atoms_has_mass(func):
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
    def wrapper(*args, **kwargs):

        if not all([atom.mass is not None for atom in args[0].atoms]):
            raise AtomicSystemMassError()

        return func(*args, **kwargs)

    return wrapper
