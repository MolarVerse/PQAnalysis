# def check_atoms_pos(func):
#     """
#     Decorator which checks that the number of atoms is equal to the number of positions.

#     Parameters
#     ----------
#     func : function
#         Function to be decorated.

#     Raises
#     ------
#     ValueError
#         If the number of atoms is not equal to the number of positions.
#     """
#     def wrapper(*args, **kwargs):
#         if args[0].pos.shape[0] != args[0].n_atoms:
#             raise ValueError(
#                 "AtomicSystem contains a different number of atoms to positions.")
#         return func(*args, **kwargs)
#     return wrapper


# def check_atoms_has_mass(func):
#     """
#     Decorator which checks that all atoms have mass information.

#     Parameters
#     ----------
#     func : function
#         Function to be decorated.

#     Raises
#     ------
#     ValueError
#         If any atom does not have mass information.
#     """
#     def wrapper(*args, **kwargs):
#         if not all([atom.mass is not None for atom in args[0].atoms]):
#             raise ValueError(
#                 "AtomicSystem contains atoms without mass information.")
#         return func(*args, **kwargs)
#     return wrapper
