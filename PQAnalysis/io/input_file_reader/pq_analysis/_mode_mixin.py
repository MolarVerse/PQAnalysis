"""
A module containing a Mixin class for modi 
related keywords of a PQAnalysis input file.
"""


from ._parse import _parse_string


class _ModeMixin:

    """
    A mixin class to read all file related
    keywords from the input dictionary.

    The following keywords are read:
        - mode
    """

    @property
    def mode(self) -> str | None:
        """str | None: Define the analysis mode."""
        return _parse_string(self.dictionary, self.mode_key)
