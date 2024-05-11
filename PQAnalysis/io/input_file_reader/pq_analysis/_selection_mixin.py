"""
A module containing a mixin class to read all selection 
related keywords from the input dictionary
"""

from ._parse import _parse_bool, _parse_string



class _SelectionMixin:

    """
    A mixin class to read all selection related keywords from the input dictionary.

    The following keywords are read:
        - selection
        - reference_selection
        - target_selection
        - use_full_atom_info
    """

    @property
    def selection(self) -> str | None:
        """str | None: The selection of the simulation."""
        return _parse_string(self.dictionary, self.selection_key)

    @property
    def reference_selection(self) -> str | None:
        """str | None: The reference selection of the simulation."""
        return _parse_string(self.dictionary, self.reference_selection_key)

    @property
    def target_selection(self) -> str | None:
        """str | None: The target selection of the simulation."""
        return _parse_string(self.dictionary, self.target_selection_key)

    @property
    def use_full_atom_info(self) -> bool | None:
        """bool | None: Whether the full atom information should be used."""
        return _parse_bool(self.dictionary, self.use_full_atom_info_key)

    @property
    def no_intra_molecular(self) -> bool | None:
        """bool | None: Whether intra molecular interactions should be excluded."""
        return _parse_bool(self.dictionary, self.no_intra_molecular_key)
