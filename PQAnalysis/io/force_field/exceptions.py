"""
Exceptions of the force-field import (:py:func:`parm2pq`).
"""

from PQAnalysis.exceptions import PQException



class Parm2PQError(PQException):

    """
    Raised if a force field cannot be converted into PQ setup files.
    """



class UnsupportedForceFieldTermError(Parm2PQError):

    """
    Raised if a force field contains a term that PQ cannot evaluate.
    """
