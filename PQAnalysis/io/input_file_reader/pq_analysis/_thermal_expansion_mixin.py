"""
A module containing a Mixin class for finite difference
related keywords of a PQAnalysis input file.
"""

from beartype.typing import List, Tuple

from ._parse import _parse_string, _parse_lists


class _ThermalExpansionMixin:
    """
    A mixin class to read all thermal expansion coefficient analysis related keywords from the input dictionary.

    The following keywords are read:
    - temperature_points
    - unit
    """

    @property
    def temperature_points(self) -> List[float] | None:
        """
        List[float] | None: The temperature points.
        """
        return _parse_lists(self.dictionary, self.temperature_points_key)

    @property
    def unit(self) -> str | None:
        """
        str | None: The unit of the box files.
        """
        return _parse_string(self.dictionary, self.unit_key)
