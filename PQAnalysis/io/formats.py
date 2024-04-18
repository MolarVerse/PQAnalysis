"""
A module containing different formats related to the io subpackage.
"""

from beartype.typing import Any

from .exceptions import BoxFileFormatError, FileWritingModeError, OutputFileFormatError
from PQAnalysis.formats import BaseEnumFormat


class OutputFileFormat(BaseEnumFormat):
    """
    An enumeration of the supported output file formats.
    """
    file_extensions = {}
    #: inference of the file format from the file extension
    AUTO = "auto"

    #: The xyz file format.
    XYZ = "xyz"
    file_extensions[XYZ] = [".xyz", ".coord", ".coords"]

    #: The vel file format.
    VEL = "vel"
    file_extensions[VEL] = [".vel", ".velocs", ".velocity"]

    #: The force file format.
    FORCE = "force"
    file_extensions[FORCE] = [".force", ".frc", ".forces"]

    #: The charge file format.
    CHARGE = "charge"
    file_extensions[CHARGE] = [".charge", ".chrg", ".charges"]

    #: The restart file format.
    RESTART = "restart"
    file_extensions[RESTART] = [".rst", ".restart"]

    #: The ENERGY file format.
    ENERGY = "energy"
    file_extensions[ENERGY] = [".en", ".energy", ".energies"]

    #: The INSTANTANEOUS ENERGY file format.
    INSTANTANEOUS_ENERGY = "instant_en"
    file_extensions[INSTANTANEOUS_ENERGY] = [
        ".instant_en", ".instant_energies", ".inst_energy"]

    #: The STRESS file format.
    STRESS = "stress"
    file_extensions[STRESS] = [".stress", ".stresses"]

    #: The VIRIAL file format.
    VIRIAL = "virial"
    file_extensions[VIRIAL] = [".virial", ".virials", ".vir"]

    #: The Info file format.
    INFO = "info"
    file_extensions[INFO] = [".info", ".information"]

    @classmethod
    def _missing_(cls, values: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        values : Tuple[Any | str] | Any
            The value to be converted to the enumeration or a tuple containing the value and the filename.

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
            raise OutputFileFormatError(
                "The file format could not be inferred from the file extension because no filename was given."
            )

        if output_file_format == cls.AUTO:
            return cls.infer_format_from_extension(filename)

        return output_file_format

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

        file_extension = file_path.split(".")[-1]

        if file_extension in cls.file_extensions[cls.XYZ]:
            return cls.XYZ
        elif file_extension in cls.file_extensions[cls.VEL]:
            return cls.VEL
        elif file_extension in cls.file_extensions[cls.FORCE]:
            return cls.FORCE
        elif file_extension in cls.file_extensions[cls.CHARGE]:
            return cls.CHARGE
        elif file_extension in cls.file_extensions[cls.RESTART]:
            return cls.RESTART
        else:
            raise OutputFileFormatError(
                f"Could not infer the file format from the file extension of \"{
                    file_path}\". Possible file formats are: {cls._member_names_}"
            )

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
