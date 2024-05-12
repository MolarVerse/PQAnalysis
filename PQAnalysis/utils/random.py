"""
A module containing functions related to generate random numbers.
"""

import secrets
import os



def get_random_seed() -> int:
    """
    Get a random seed for the simulation.

    Returns
    -------
    int:
        The random seed. If the environment variable CONSTANT_SEED is set, 
        the seed will be the value of the environment variable. Otherwise,
        the seed will be a random 128-bit integer.
    """
    if 'CONSTANT_SEED' in os.environ:
        return int(os.environ['CONSTANT_SEED'])

    return secrets.randbits(128)
