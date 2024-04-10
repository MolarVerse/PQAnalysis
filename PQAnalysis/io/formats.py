"""
A module containing different formats related to the io subpackage.
"""

from beartype.typing import Any

from .exceptions import BoxFileFormatError, FileWritingModeError, OutputFileFormatError
from PQAnalysis.formats import BaseEnumFormat


class OutputFileFormat(BaseEnumFormat):
    #: inference of the file format from the file extension
    AUTO = "auto"

    #: The xyz file format.
    XYZ = "xyz"

    #: The restart file format.
    RESTART = "restart"

    @classmethod
    def _missing_(cls, value: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        value : Any
            The value to return.

        Returns
        -------
        Any
            The value to return.
        """

        return super()._missing_(value, OutputFileFormatError)

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

        if file_path.endswith(".xyz"):
            return cls.XYZ
        elif file_path.endswith(".rst"):
            return cls.RESTART
        else:
            raise OutputFileFormatError(
                f"Could not infer the file format from the file extension of \"{file_path}\". Possible file formats are: {cls._member_names_}"
            )


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

    @ classmethod
    def _missing_(cls, value: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        value : Any
            The value to return.

        Returns
        -------
        Any
            The value to return.
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
    #: | where all X represent the vertices of the box. The first line contains the number of vertices.
    #: | The second line contains the box dimensions and box angles as the comment line for a xyz file.
    VMD = "vmd"

    #: | The data file format.
    #: | The format looks in general like this:
    #: |            1 1.0 1.0 1.0 90.0 90.0 90.0
    #: |            2 1.0 1.0 1.0 90.0 90.0 90.0
    #: |            ...
    #: |            n 1.1 1.1 1.1 90.0 90.0 90.0
    #: | where the first column represents the step starting from 1, the second to fourth column
    #: | represent the box vectors a, b, c, the fifth to seventh column represent the box angles.
    DATA = "data"

    @ classmethod
    def _missing_(cls, value: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        value : Any
            The value to return.

        Returns
        -------
        Any
            The value to return.
        """

        return super()._missing_(value, BoxFileFormatError)
