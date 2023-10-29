"""
A module containing the base class for all writers.

...

Classes
-------
BaseWriter
    A base class for all writers.
"""

import sys
import os


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

    Methods
    -------
    open()
        Opens the file to write to.
    close()
        Closes the file to write to.
    """

    def open(self):
        """
        Opens the file to write to.
        """
        if self.filename is not None:
            self.file = open(self.filename, self.mode)

    def close(self):
        """
        Closes the file to write to.
        """
        if self.file is not None:
            self.file.close()

    def __init__(self, filename: str, mode: str):
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
            raise ValueError('Invalid mode has to be either \'w\' or \'a\'.')
        elif mode == 'w' and filename is not None and os.path.isfile(filename):
            raise ValueError(
                f"File {filename} already exists. Use mode \'a\' to append to file.")

        if filename is None:
            self.file = sys.stdout
        else:
            self.file = None

        self.mode = 'a'
        self.filename = filename
        self.format = format
