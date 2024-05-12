"""
A module containing a Mixin class for position
related keywords of a PQAnalysis input file.
"""

from PQAnalysis.types import PositiveInt, PositiveReal
from ._parse import _parse_positive_int, _parse_positive_real



class _PositionsMixin:

    """
    A mixin class to read all position related keywords from the input dictionary.

    The following keywords are read:
        - r_max
        - r_min
        - delta_r
        - n_bins
    """

    @property
    def r_max(self) -> PositiveReal | None:
        """PositiveReal | None: The maximum radius of the PQAnalysis."""
        return _parse_positive_real(self.dictionary, self.r_max_key)

    @property
    def r_min(self) -> PositiveReal | None:
        """PositiveReal | None: The minimum radius of the PQAnalysis."""
        return _parse_positive_real(self.dictionary, self.r_min_key)

    @property
    def delta_r(self) -> PositiveReal | None:
        """PositiveReal | None: The radius step size of the PQAnalysis."""
        return _parse_positive_real(self.dictionary, self.delta_r_key)

    @property
    def n_bins(self) -> PositiveInt | None:
        """PositiveInt | None: The number of bins of the PQAnalysis."""
        return _parse_positive_int(self.dictionary, self.n_bins_key)
