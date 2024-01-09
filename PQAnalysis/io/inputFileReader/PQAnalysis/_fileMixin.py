from beartype.typing import List

from ._parse import _parse_files, _parse_string


class _FileMixin:
    """
    A mixin class to read all file related keywords from the input dictionary.

    The following keywords are read:
        - traj_files
        - out_file
        - log_file
    """
    @property
    def traj_files(self) -> List[str] | None:
        """
        Returns the trajectory files from the input dictionary.

        Returns
        -------
        List[str] | None
            the trajectory files or None if the key is not in the dictionary
        """
        return _parse_files(self.dictionary, self.traj_files_key)

    @property
    def out_file(self) -> str | None:
        """
        Returns the output file from the input dictionary.

        Returns
        -------
        str | None
            the output file or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.out_file_key)

    @property
    def log_file(self) -> str | None:
        """
        Returns the log file from the input dictionary.

        Returns
        -------
        str | None
            the log file or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.log_file_key)

    @property
    def moldescriptor_file(self) -> str | None:
        """
        Returns the moldescriptor file from the input dictionary.

        Returns
        -------
        str | None
            the moldescriptor file or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.moldescriptor_file_key)

    @property
    def restart_file(self) -> str | None:
        """
        Returns the restart file from the input dictionary.

        Returns
        -------
        str | None
            the restart file or None if the key is not in the dictionary
        """
        return _parse_string(self.dictionary, self.restart_file_key)
