"""
A module containing different formats related to the io subpackage.
"""
import logging

from beartype.typing import Any, List

from PQAnalysis.formats import BaseEnumFormat
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from .exceptions import (
    BoxFileFormatError, FileWritingModeError, OutputFileFormatError
)

logger = logging.getLogger(__package_name__).getChild("OutputFileFormat")
logger = setup_logger(logger)



class OutputFileFormat(BaseEnumFormat):

    """
    An enumeration of the supported output file formats.
    """
    #: inference of the file format from the file extension
    AUTO = "auto"

    #: The xyz file format.
    XYZ = "xyz"

    #: The vel file format.
    VEL = "vel"

    #: The force file format.
    FORCE = "force"

    #: The charge file format.
    CHARGE = "charge"

    #: The restart file format.
    RESTART = "restart"

    #: The ENERGY file format.
    ENERGY = "energy"

    #: The INSTANTANEOUS ENERGY file format.
    INSTANTANEOUS_ENERGY = "instant_en"

    #: The STRESS file format.
    STRESS = "stress"

    #: The VIRIAL file format.
    VIRIAL = "virial"

    #: The Info file format.
    INFO = "info"

    @classmethod
    def _missing_(cls, values: Any) -> Any:  # pylint: disable=arguments-differ
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        values : Tuple[Any | str] | Any
            The value to be converted to the enumeration
            or a tuple containing the value and the filename.

        Returns
        -------
        Any
            The value to return.
        """

        if isinstance(values, tuple):
            value = values[0]
            filename = values[1]
        else:
            value = values
            filename = None

        output_file_format = super()._missing_(value, OutputFileFormatError)

        if output_file_format == cls.AUTO and filename is None:
            logger.error(
                (
                    "The file format could not be inferred from "
                    "the file extension because no filename was given."
                ),
                exception=OutputFileFormatError
            )

        if output_file_format == cls.AUTO:
            return cls.infer_format_from_extension(filename)

        return output_file_format

    @classmethod
    def file_extensions(cls):
        """
        Get the file extensions of the different file formats.

        Returns
        -------
        dict[str, list[str]]
            The file extensions of the different file formats.
        """
        file_extensions = {}
        file_extensions[cls.XYZ.value] = [".xyz", ".coord", ".coords"]
        file_extensions[cls.VEL.value] = [".vel", ".velocs", ".velocity"]
        file_extensions[cls.FORCE.value] = [".force", ".frc", ".forces"]
        file_extensions[cls.CHARGE.value] = [".charge", ".chrg", ".charges"]
        file_extensions[cls.RESTART.value] = [".rst", ".restart"]
        file_extensions[cls.ENERGY.value] = [".en", ".energy", ".energies"]
        file_extensions[cls.STRESS.value] = [".stress", ".stresses"]
        file_extensions[cls.VIRIAL.value] = [".virial", ".virials", ".vir"]
        file_extensions[cls.INFO.value] = [".info", ".information"]
        file_extensions[cls.INSTANTANEOUS_ENERGY.value] = [
            ".instant_en", ".instant_energies", ".inst_energy"
        ]

        return file_extensions

    @classmethod
    def infer_format_from_extension(cls, file_path: str) -> "OutputFileFormat":
        """
        Infer the file format from the file extension.

        Parameters
        ----------
        file_path : str
            The file path to infer the file format from.

        Returns
        -------
        OutputFileFormat
            The inferred file format.
        """

        file_extension = "." + file_path.split(".")[-1]

        if file_extension in cls.file_extensions()[cls.XYZ.value]:
            return cls.XYZ

        if file_extension in cls.file_extensions()[cls.VEL.value]:
            return cls.VEL

        if file_extension in cls.file_extensions()[cls.FORCE.value]:
            return cls.FORCE

        if file_extension in cls.file_extensions()[cls.CHARGE.value]:
            return cls.CHARGE

        if file_extension in cls.file_extensions()[cls.RESTART.value]:
            return cls.RESTART

        logger.error(
            (
                "Could not infer the file format from the file extension of "
                f"\"{file_path}\". Possible file formats are: "
                f"{cls.__members__.values()}"
            ),
            exception=OutputFileFormatError
        )

    @classmethod
    def get_file_extensions(
        cls, file_format: "OutputFileFormat | str"
    ) -> List[str]:
        """
        Get the file extensions of the given file format.

        Parameters
        ----------
        file_format : OutputFileFormat
            The file format to get the file extensions for.

        Returns
        -------
        list[str]
            The file extensions of the given file format.
        """

        file_format = cls(file_format)

        return cls.file_extensions()[file_format.value]

    @classmethod
    def find_matching_files(
        cls,
        file_path: List[str],
        output_file_format: "OutputFileFormat | str",
        extension: str | None = None
    ) -> List[str]:
        """
        Find the files that match the given file format.

        Parameters
        ----------
        file_path : List[str]
            The file paths to search for the files.
        output_file_format : OutputFileFormat | str
            The file format to search for.
        extension : str | None, optional
            The extension to search for, by default None.
            If None, all files with the given file format are returned.
            Else, only the files with the given extension are returned.

        Returns
        -------
        List[str]
            The files that match the given file format.
        """
        if extension is not None:
            files = [file for file in file_path if file.endswith(extension)]
        else:
            files = [
                file for file in file_path if file.
                endswith(tuple(cls.get_file_extensions(output_file_format)))
            ]

        return files

    def lower(self) -> str:
        """
        This method returns the lower case representation of the enumeration.

        Returns
        -------
        str
            The lower case representation of the enumeration.
        """

        return self.value.lower()

    def __eq__(self, other: object) -> bool:
        """
        This method checks if the enumeration is equal to the other object.

        Parameters
        ----------
        other : object
            The object to compare with.

        Returns
        -------
        bool
            True if the enumeration is equal to the other object, False otherwise.
        """
        return super().__eq__(other) or self.value == str(other)



class FileWritingMode(BaseEnumFormat):

    """
    An enumeration of the supported file write modes.

    """

    #: The write mode for overwriting a file.
    OVERWRITE = "o"

    #: The write mode for appending to a file.
    APPEND = "a"

    #: The write mode for writing to a file
    WRITE = "w"

    @classmethod
    def _missing_(cls, value: Any) -> Any:  # pylint: disable=arguments-differ
        """
        This method returns the missing value of the enumeration.

        Parameters
        - ---------
        value: Any
            The value to return .

        Returns
        - ------
        Any
            The value to return .
        """

        return super()._missing_(value, FileWritingModeError)



class BoxFileFormat(BaseEnumFormat):

    """
    An enumeration of the supported box file formats.

    """

    #: | The VMD format.
    #: | The format looks in general like this:
    #: |            8
    #: |            Box  1.0  1.0  1.0    90.0 90.0 90.0
    #: |            X   -0.5 -0.5 -0.5
    #: |            X    0.5 -0.5 -0.5
    #: |            X   -0.5  0.5 -0.5
    #: |            X    0.5  0.5 -0.5
    #: |            X   -0.5 -0.5  0.5
    #: |            X    0.5 -0.5  0.5
    #: |            X   -0.5  0.5  0.5
    #: |            X    0.5  0.5  0.5
    #: |            8
    #: |            Box  2.0  2.0  2.0    90.0 90.0 90.0
    #: |            X   -1.0 -1.0 -1.0
    #: |            ...
    #: | where all X represent the vertices of the box.
    # The first line contains the number of vertices.
    #: | The second line contains the box dimensions and
    # box angles as the comment line for a xyz file.
    VMD = "vmd"

    #: | The data file format.
    #: | The format looks in general like this:
    #: |            1 1.0 1.0 1.0 90.0 90.0 90.0
    #: |            2 1.0 1.0 1.0 90.0 90.0 90.0
    #: |            ...
    #: |            n 1.1 1.1 1.1 90.0 90.0 90.0
    #: | where the first column represents the step
    # starting from 1, the second to fourth column
    #: | represent the box vectors a, b, c,
    # the fifth to seventh column represent the box angles.
    DATA = "data"

    @classmethod
    def _missing_(cls, value: Any) -> Any:  # pylint: disable=arguments-differ
        """
        This method returns the missing value of the enumeration.

        Parameters
        - ---------
        value: Any
            The value to return .

        Returns
        - ------
        Any
            The value to return .
        """

        return super()._missing_(value, BoxFileFormatError)
