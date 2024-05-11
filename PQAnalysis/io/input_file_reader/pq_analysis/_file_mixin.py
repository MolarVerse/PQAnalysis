"""
A module containing a Mixin class for file 
related keywords of a PQAnalysis input file.
"""

from beartype.typing import List

from ._parse import _parse_files, _parse_string



class _FileMixin:

    """
    A mixin class to read all file related
    keywords from the input dictionary.

    The following keywords are read:
        - traj_files
        - out_file
        - log_file
    """

    @property
    def traj_files(self) -> List[str] | None:
        """List[str] | None: The trajectory files of the simulation."""
        return _parse_files(self.dictionary, self.traj_files_key)

    @property
    def out_file(self) -> str | None:
        """str | None: The out file of the simulation."""
        return _parse_string(self.dictionary, self.out_file_key)

    @property
    def log_file(self) -> str | None:
        """str | None: The log file of the simulation."""
        return _parse_string(self.dictionary, self.log_file_key)

    @property
    def moldescriptor_file(self) -> str | None:
        """str | None: The moldescriptor file of the simulation."""
        return _parse_string(self.dictionary, self.moldescriptor_file_key)

    @property
    def restart_file(self) -> str | None:
        """str | None: The restart file of the simulation."""
        return _parse_string(self.dictionary, self.restart_file_key)
