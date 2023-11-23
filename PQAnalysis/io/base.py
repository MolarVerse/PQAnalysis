"""
A module containing the base class for all writers and readers.

...

Classes
-------
BaseWriter
    A base class for all writers.
BaseReader
    A base class for all readers.
"""

import sys
import os

from beartype.typing import List


class BaseWriter:
    """
    A base class for all writers.

    ...

    Attributes
    ----------
    file : file
        The file to write to.
    mode : str
        The mode of the file. Either 'w' for write or 'a' for append.
    filename : str
        The name of the file to write to. If None, the output is printed to stdout.
    """

    def __init__(self, filename: str | None = None, mode: str | None = 'w') -> None:
        """
        It sets the file to write to - either a file or stdout (if filename is None) - and the mode of the file.

        Parameters
        ----------
        filename : str
            The name of the file to write to. If None, the output is printed to stdout.
        mode : str
            The mode of the file. Either 'w' for write or 'a' for append.

        Raises
        ------
        ValueError
            If the given mode is not 'w' or 'a'.
        ValueError
            If the given filename already exists and the mode is 'w'.
        """
        if mode not in ['w', 'a']:
            raise ValueError('Invalid mode - has to be either \'w\' or \'a\'.')
        elif mode == 'w' and filename is not None and os.path.isfile(filename):
            raise ValueError(
                f"File {filename} already exists. Use mode \'a\' to append to file.")

        if filename is None:
            self.file = sys.stdout
        else:
            self.file = None

        self.mode = 'a'
        self.filename = filename

    def open(self) -> None:
        """
        Opens the file to write to.
        """
        if self.filename is not None:
            self.file = open(self.filename, self.mode)

    def close(self) -> None:
        """
        Closes the file to write to.
        """
        if self.file is not None and self.filename is not None:
            self.file.close()

        self.file = None


class BaseReader:
    """
    A base class for all readers.

    ...

    Attributes
    ----------
    filename : str
        The name of the file to read from.
    filenames : list of str
        The filenames to read from.
    multiple_files : bool
        Whether the reader reads from a single file or multiple files.
    """

    def __init__(self, filename: str | List[str]) -> None:
        """
        Initializes the BaseReader with the given filename.

        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames to read from.

        Raises
        ------
        FileNotFoundError
            If the given file does not exist.
        """

        if isinstance(filename, str):
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"File {filename} not found.")

            self.filename = filename
            self.multiple_files = False
        else:
            if not all([os.path.isfile(f) for f in filename]):
                raise FileNotFoundError(
                    f"At least one of the given files does not exist.")

            self.filenames = filename
            self.multiple_files = True
