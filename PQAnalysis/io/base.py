"""
A module containing the base class for all writers and readers.
"""

import sys
import os

from beartype.typing import List


class BaseWriter:
    """
    A base class for all writers.

    This class can be used as a base class for all writers. It provides the functionality to write to a file or stdout. Meaning, if the filename is None, the output is printed to stdout. Otherwise, the output is written to the given file. The mode of the file can be either 'w' for write or 'a' for append. If the mode is 'w' and the file already exists, a ValueError is raised. If the mode is 'a' and the file does not exist, it is created. If the mode is 'o' and the file already exists, it is overwritten.

    Attention
    ---------

    Please be aware that after the initialization of the BaseWriter, the file is not opened yet. To open the file, the open() method has to be called. To close the file, the close() method has to be called. After calling the open() method of a file that already exists, the mode is automatically set to 'a'. This means that for all subsequent calls of the open() method, the mode is 'a' and not 'w'. This is to prevent overwriting the file by accident. 

    """

    #: The modes of the file.
    #: 'w' for write, 'a' for append and 'o' for overwrite.
    modes = ['w', 'a', 'o']

    def __init__(self, filename: str | None = None, mode: str | None = 'w') -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the file to write to. If None, the output is printed to stdout.
        mode : str
            The mode of the file. Either 'w' for write or 'a' for append.

        Raises
        ------
        ValueError
            If the given mode is not in .
        ValueError
            If the given filename already exists and the mode is 'w'.
        """

        if filename is None:
            self.file = sys.stdout
        else:
            self.file = None

        self.filename = filename

        self.mode = mode

    def open(self) -> None:
        """
        Opens the file to write to.
        """
        if self.filename is not None:
            self.file = open(self.filename, self.mode)
            self.mode = 'a'

    def close(self) -> None:
        """
        Closes the file to write to.
        """
        if self.file is not None and self.filename is not None:
            self.file.close()

        self.file = None

    @property
    def mode(self) -> str:
        """
        str: The mode of the file.

        The setter checks if the given mode is in BaseWriter.modes. Furthermore, if the mode is 'w' and the file already exists, a ValueError is raised. If the mode is 'a' and the file does not exist, it is created. If the mode is 'o' and the file already exists, it is overwritten.
        """
        return self._mode

    @mode.setter
    def mode(self, mode: str) -> None:
        if mode not in self.modes:
            raise ValueError(
                'Invalid mode - has to be one of ' + str(self.modes))

        if mode not in self.modes:
            raise ValueError(
                'Invalid mode - has to be one of ' + str(self.modes))

        if mode == 'w' and self.filename is not None and os.path.isfile(self.filename):
            raise ValueError(
                f"File {self.filename} already exists. Use mode \'a\' to append to the file or mode \'o\' to overwrite the file.")

        self._mode = mode


class BaseReader:
    """
    A base class for all readers.

    This class can be used as a base class for all Readers implemented in the PQAnalysis package. It provides the functionality to read from a file or multiple files. If the filename is a string, the file is read from. The basic idea is just to ease the handling of multiple files and to provide a common interface for all readers. The actual reading is done in the methods of the subclasses. This wrapper only sets the filename/filenames and checks if the file/files exist. Furthermore, a boolean is set to indicate if the reader reads from a single file or multiple files.
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
            filenames = filename
            for filename in filenames:
                if not os.path.isfile(filename):
                    raise FileNotFoundError(
                        f"At least one of the given files does not exist. File {filename} not found.")

            self.filenames = filenames
            self.multiple_files = True
