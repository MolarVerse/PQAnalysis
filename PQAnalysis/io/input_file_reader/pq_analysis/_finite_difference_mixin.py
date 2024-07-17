"""
A module containing a Mixin class for finite difference
related keywords of a PQAnalysis input file.
"""

from beartype.typing import List, Tuple

from ._parse import _parse_lists


class _FiniteDifferenceMixin:
    """
    A mixin class to read all finite difference related keywords from the input dictionary.

    The following keywords are read:
    - finite_difference_points
    - std_points
    - temperature_points
    """

    @property
    def finite_difference_points(self) -> Tuple[List[float]] | List[float] | None:
        """
        Tuple[List[float]] | List[float]| None: The finite difference points.
        It is a list of floats or a tuple of lists of floats. 

        """
        return _parse_lists(self, "finite_difference_points")

    @property
    def temperature_points(self) -> List[float] | None:
        """
        List[float] | None: The temperature points.
        """
        return _parse_lists(self, "temperature_points")

    @property
    def std_points(self) -> List[float] | None:
        """
        List[float] | None: The standard deviation points.
        """
        return _parse_lists(self, "std_points")
