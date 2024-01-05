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
        """
        Returns the selection from the input dictionary as a string for the Lark selection parser.

        Returns
        -------
        str | None
            the selection or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.selection_key)

    @property
    def reference_selection(self) -> str | None:
        """
        Returns the reference selection from the input dictionary as a string for the Lark selection parser.

        Returns
        -------
        str | None
            the reference selection or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.reference_selection_key)

    @property
    def target_selection(self) -> str | None:
        """
        Returns the target selection from the input dictionary as a string for the Lark selection parser.

        Returns
        -------
        str | None
            the target selection or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.target_selection_key)

    @property
    def use_full_atom_info_for_selection(self) -> bool | None:
        """
        Returns the use_full_atom_info from the input dictionary as a bool.

        This information is used to determine if the full atom information should be used for the selection or not.

        Returns
        -------
        bool | None
            the use_full_atom_info or None if the key is not in the dictionary
        """
        return _parse_bool(self.dictionary, self.use_full_atom_info_key)
