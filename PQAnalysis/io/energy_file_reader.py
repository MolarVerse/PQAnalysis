"""
A module containing the EnergyFileReader class.
"""

import os
import logging
import numpy as np

from PQAnalysis.physical_data import Energy
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQFileNotFoundError
from PQAnalysis.type_checking import runtime_type_checking

from .base import BaseReader
from .info_file_reader import InfoFileReader



class EnergyFileReader(BaseReader):

    """
    A class to read energy files from molecular dynamics simulations.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str,
        info_filename: str | None = None,
        use_info_file: bool = True,
        engine_format: MDEngineFormat | str = MDEngineFormat.PQ
    ) -> None:
        """
        For the initialization of the EnergyFileReader, the filename of the energy
        file has to be given. The info_filename can be given as well. If no
        info_filename is given, the energy filename is used to find the info 
        file with the proper extension and the base name of the energy file. 
        If a info_filename is given, this filename is used to find the info 
        file. If the info_filename was explicitly set to a non-existing file,
        a FileNotFoundError is raised. If use_info_file is set to False, no info
        file is searched for.

        Providing info files is optional, but where possible, it is recommended to 
        provide an info file. The info file contains information about the physical
        quantities stored in the energy file. If an info file is provided, the 
        physical quantities are assumed to be in the order of the info file. The 
        info file can also contain units for the physical quantities. If units are
        provided, the units are stored in the Energy object, which is returned by 
        the read method.

        For more information about the energy object, see
        :py:class:`~PQAnalysis.physicalData.energy.Energy`.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        info_filename : str, optional
            The name of the info file to read from, by default None
        use_info_file : bool, optional
            If True, the info file is searched for, by default True
        engine_format : MDEngineFormat | str, optional
            The format of the file, by default MDEngineFormat.PQ
        """
        super().__init__(filename)
        self.info_filename = info_filename

        if use_info_file:
            self.with_info_file = self.__info_file_found__()
        else:
            self.with_info_file = False

        self.format = MDEngineFormat(engine_format)

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

        if self.with_info_file:

            reader = InfoFileReader(
                self.info_filename,
                engine_format=self.format
            )

            info, units = reader.read()

        with open(self.filename, "r", encoding='utf-8') as file:

            data = []

            for line in file:
                if line.startswith("#"):
                    continue

                data_line = [float(x) for x in line.split()]
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
            except PQFileNotFoundError:
                self.info_filename = None

        else:

            try:
                BaseReader(self.info_filename)
            except PQFileNotFoundError:
                self.logger.error(
                    f"Info File {self.info_filename} not found.",
                    exception=PQFileNotFoundError
                )

        return self.info_filename is not None
