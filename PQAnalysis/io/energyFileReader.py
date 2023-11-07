"""
A module containing the EnergyFileReader class.

...

Classes
-------
EnergyFileReader
    A class to read energy files.
"""

import os
import numpy as np

from .base import BaseReader
from .infoFileReader import InfoFileReader
from ..physicalData.energy import Energy


class EnergyFileReader(BaseReader):
    """
    A class to read energy files.

    Parameters
    ----------
    BaseReader : BaseReader
        A base class for all readers.

    Attributes
    ----------
    filename : str
        The name of the file to read from.
    info_filename : str
        The name of the info file to read from.
    withInfoFile : bool
        If True, the info file was found.
    """

    formats = ["pimd-qmcf", "qmcfc"]

    def __init__(self,
                 filename: str,
                 info_filename: str | None = None,
                 use_info_file: bool = True,
                 format: str = "pimd-qmcf"
                 ) -> None:
        """
        Initializes the EnergyFileReader with the given filename.

        If no info_filename is given, the energy filename is used to find the
        info file. If a info_filename is given, this filename is used to find the
        info file. If the info_filename was explicitly set to a non-existing file,
        a FileNotFoundError is raised. If use_info_file is set to False, no info
        file is searched for.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        info_filename : str, optional
            The name of the info file to read from, by default None
        use_info_file : bool, optional
            If True, the info file is searched for, by default True
        format : str, optional
            The format of the info file, by default "pimd-qmcf"

        Raises
        ------
        ValueError
            If the format is not supported.
        """
        super().__init__(filename)
        self.info_filename = info_filename

        if use_info_file:
            self.withInfoFile = self.__info_file_found__()
        else:
            self.withInfoFile = False

        if format not in self.formats:
            raise ValueError(
                f"Format {format} is not supported. Supported formats are {self.formats}.")

        self.format = format

    def read(self) -> Energy:
        """
        Reads the energy file.

        The data of the energy file is returned as a np.array within a Energy object.
        The data array is stored in a way that each column corresponds to a physical
        quantity. The order of the columns is the same as in the info file.

        Returns
        -------
        Energy
            The data of the energy file as a np.array within a Energy object.
            In addition, the info and units of the info file are stored in the Energy
            object, if an info file was found.
        """
        info, units = None, None

        if self.withInfoFile:
            reader = InfoFileReader(self.info_filename, format=self.format)
            info, units = reader.read()

        with open(self.filename, "r") as file:

            data = []

            for line in file:
                if line.startswith("#"):
                    continue

                data_line = map(lambda x: float(x), line.split())
                data.append(list(data_line))

        return Energy(np.array(data).T, info, units)

    def __info_file_found__(self) -> bool:
        """
        Checks if a info file exists for the given file.

        If no info_filename is given, the energy filename is used to find the
        info file. If a info_filename is given, this filename is used to find the
        info file. If the info_filename was explicitly set to a non-existing file,
        a FileNotFoundError is raised.

        Returns
        -------
        bool
            True if a info file was found, False otherwise.

        Raises
        ------
        FileNotFoundError
            If an explicitly given info file does not exist.
        """
        if self.info_filename is None:

            self.info_filename = os.path.splitext(self.filename)[0] + ".info"
            try:
                BaseReader(self.info_filename)
            except FileNotFoundError:
                self.info_filename = None
        else:
            try:
                BaseReader(self.info_filename)
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Info File {self.info_filename} not found.")

        if self.info_filename is not None:
            print(f"A Info File \'{self.info_filename}\' was found.")
            return True
        else:
            return False
