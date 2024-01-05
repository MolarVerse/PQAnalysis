from ._parse import _parse_positive_int, _parse_positive_real
from ....types import PositiveInt, PositiveReal


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
        """
        Returns the r_max from the input dictionary as a positive real number.

        Returns
        -------
        PositiveReal | None
            the r_max or None if the key is not in the dictionary
        """
        return _parse_positive_real(self.dictionary, self.r_max_key)

    @property
    def r_min(self) -> PositiveReal | None:
        """
        Returns the r_min from the input dictionary as a positive real number.

        Returns
        -------
        PositiveReal | None
            the r_min or None if the key is not in the dictionary
        """
        return _parse_positive_real(self.dictionary, self.r_min_key)

    @property
    def delta_r(self) -> PositiveReal | None:
        """
        Returns the delta_r from the input dictionary as a positive real number.

        Returns
        -------
        PositiveReal | None
            the delta_r or None if the key is not in the dictionary
        """
        return _parse_positive_real(self.dictionary, self.delta_r_key)

    @property
    def n_bins(self) -> PositiveInt | None:
        """
        Returns the number of bins from the input dictionary as a positive integer.

        Returns
        -------
        PositiveInt | None
            the number of bins or None if the key is not in the dictionary
        """
        return _parse_positive_int(self.dictionary, self.n_bins_key)
