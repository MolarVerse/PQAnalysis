"""
A module containing the EnergyFileReader class.
"""

import os
import logging
import math
import re
from numbers import Real
import numpy as np

from PQAnalysis.physical_data import Energy
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQFileNotFoundError
from PQAnalysis.type_checking import runtime_type_checking

from .base import BaseReader
from .energy_file_header import EnergyLabels, read_energy_header
from .exceptions import EnergyFileReaderError
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
        filename: str | list[str],
        info_filename: str | None = None,
        use_info_file: bool = True,
        engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
        timestep_fs: Real | None = None
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
        filename : str or list of str
            The name of the file to read from or a list of filenames to read
            from.
        info_filename : str, optional
            The name of the info file to read from, by default None
        use_info_file : bool, optional
            If True, the info file is searched for, by default True.
            Without info labels, numeric columns are returned unchanged.
        engine_format : MDEngineFormat | str, optional
            The format of the file, by default MDEngineFormat.PQ
        timestep_fs : Real, optional
            Constant PQ integration timestep in femtoseconds. It must be
            finite and positive and agree with any energy header. Native
            ``# timestep = <value> fs`` headers are used automatically.
            With PQ info labels but no timestep, the first column remains
            integration steps and its unit is ``step``, not ``ps``.
            For files with the PQ column header it is only used if the
            header gives ``SIMULATION-TIME`` in ``step`` or ``STEP``
            in ``-``.

        Raises
        ------
        EnergyFileReaderError
            If the explicit timestep is invalid or used with QMCFC input.
        """
        super().__init__(filename)
        self.info_filename = info_filename

        if use_info_file:
            self.with_info_file = self.__info_file_found__()
        else:
            self.with_info_file = False

        self.format = MDEngineFormat(engine_format)
        self.timestep_fs = timestep_fs
        if timestep_fs is not None:
            if self.format != MDEngineFormat.PQ:
                raise EnergyFileReaderError(
                    "timestep_fs is only supported for PQ energy files."
                )
            self._validate_timestep(timestep_fs)

    def read(self) -> Energy:
        """
        Reads the energy file.

        The data of the energy file is returned as a np.array within a Energy object.
        The data array is stored in a way that each column corresponds to a physical
        quantity. The order of the columns is the same as in the info file.

        If the energy file starts with the PQ header
        (``# PQ energy file``), the column names and units are taken
        from its column line instead of the info file.

        Returns
        -------
        Energy
            The data of the energy file as a np.array within a Energy object.
            In addition, the info and units of the info file are stored in the Energy
            object, if an info file was found.
        """
        info, units, header = self._labels()

        pq_time_axis = (
            self.format == MDEngineFormat.PQ and info is not None and
            "SIMULATION-TIME" in info and
            (header is None or units["SIMULATION-TIME"] == "step")
        )
        data = []
        timesteps = []

        for filename in self._filenames_to_read():
            file_data = self._read_data_from_file(filename)

            if data and file_data and len(file_data[0]) != len(data[0]):
                self.logger.error(
                    (
                        f"Energy file {filename} has {len(file_data[0])} "
                        f"columns, but the previously read energy files "
                        f"have {len(data[0])} columns."
                    ),
                    exception=EnergyFileReaderError
                )

            if header is not None:
                _check_header_width(filename, file_data, header)

            timestep = self._file_timestep(filename) if pq_time_axis else None
            timesteps.append(timestep)
            if pq_time_axis and timestep is not None:
                column = info["SIMULATION-TIME"]
                for row in file_data:
                    row[column] *= timestep / 1000.0
            data.extend(file_data)

        if pq_time_axis:
            known = [value is not None for value in timesteps]
            if any(known) and not all(known):
                raise EnergyFileReaderError(
                    "Cannot combine PQ energy files with and without "
                    "timestep metadata. Provide timestep_fs or read "
                    "the files separately."
                )
            units = dict(units)
            units["SIMULATION-TIME"] = "ps" if all(known) else "step"

        return Energy(np.array(data).T, info, units)

    def _labels(self) -> tuple[dict | None, dict | None, EnergyLabels | None]:
        """Return info, units and the PQ header labels, header first."""
        info, units = None, None
        if self.with_info_file:
            info, units = InfoFileReader(
                self.info_filename, engine_format=self.format
            ).read()

        header = self._header_labels()
        if header is None:
            return info, units, None

        if info is not None and not _same_labels((info, units), header):
            self.logger.warning(
                f"The columns in the header of the energy file differ "
                f"from the info file {self.info_filename}. Using the "
                "header."
            )
        return *header, header

    def _header_labels(self) -> EnergyLabels | None:
        """Return the column labels of the PQ header, if the files have one."""
        headers = [
            header for header in (
                read_energy_header(filename)
                for filename in self._filenames_to_read()
            ) if header is not None
        ]
        if any(header != headers[0] for header in headers[1:]):
            raise EnergyFileReaderError(
                "The energy files have different column headers."
            )
        return headers[0] if headers else None

    @staticmethod
    def _validate_timestep(value: Real) -> None:
        """Reject nonphysical integration timestep metadata."""
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise EnergyFileReaderError(
                "The PQ energy timestep must be finite and positive."
            )

    def _file_timestep(self, filename: str) -> float | None:
        """Resolve one PQ file's explicit or header timestep in fs."""
        declared = None
        with open(filename, encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                if not text.startswith("#"):
                    break
                value = self._parse_timestep_header(text, filename)
                if value is None:
                    continue
                if declared is not None and not math.isclose(
                    declared, value, rel_tol=1e-12
                ):
                    raise EnergyFileReaderError(
                        f"Conflicting timestep metadata in {filename}."
                    )
                declared = value
        if self.timestep_fs is not None:
            if declared is not None and not math.isclose(
                declared, self.timestep_fs, rel_tol=1e-12
            ):
                raise EnergyFileReaderError(
                    f"timestep_fs conflicts with the header in {filename}."
                )
            return float(self.timestep_fs)
        return declared

    @classmethod
    def _parse_timestep_header(cls, text: str, filename: str) -> float | None:
        """Read a timestep declaration, ignoring unrelated comments."""
        if not re.match(r"#\s*timestep\s*=", text, re.IGNORECASE):
            return None
        match = re.fullmatch(
            r"#\s*timestep\s*=\s*(\S+)\s+(\S+)",
            text,
            re.IGNORECASE,
        )
        if match is None:
            raise EnergyFileReaderError(
                f"Invalid timestep metadata in {filename}. "
                "Expected '# timestep = <value> fs'."
            )
        try:
            value = float(match[1])
        except ValueError as error:
            raise EnergyFileReaderError(
                f"Invalid timestep metadata in {filename}."
            ) from error
        cls._validate_timestep(value)
        unit = match[2].lower()
        if unit not in ("fs", "ps"):
            raise EnergyFileReaderError(
                f"Unsupported timestep unit {match[2]} in {filename}."
            )
        value *= 1000.0 if unit == "ps" else 1.0
        cls._validate_timestep(value)
        return value

    @classmethod
    def _read_data_from_file(cls, filename: str) -> list[list[float]]:
        """
        Reads energy data from a single file.
        """
        data = []

        with open(filename, "r", encoding='utf-8') as file:
            for line_number, line in enumerate(file, start=1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue

                data_line = [float(x) for x in line.split()]

                if data and len(data_line) != len(data[0]):
                    cls.logger.error(
                        (
                            f"Invalid number of columns in energy file "
                            f"{filename} line {line_number}. Expected "
                            f"{len(data[0])} columns."
                        ),
                        exception=EnergyFileReaderError
                    )

                data.append(list(data_line))

        return data

    def _filenames_to_read(self) -> list[str]:
        """
        Gets the energy filenames to read.
        """
        if self.multiple_files:
            return self.filenames

        return [self.filename]

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

            self.info_filename = os.path.splitext(
                self._filenames_to_read()[0]
            )[0] + ".info"

            if not os.path.isfile(self.info_filename):
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



def _same_labels(info_labels, header_labels) -> bool:
    """Compare labels; the info file's time unit describes its summary."""
    info, units = info_labels
    header_info, header_units = header_labels
    if info != header_info:
        return False
    return all(
        (units or {}).get(name) == header_units[name]
        for name in info
        if name != "SIMULATION-TIME"
    )



def _check_header_width(
    filename: str,
    file_data: list[list[float]],
    header: EnergyLabels,
) -> None:
    """Raise if the data width differs from the header column count."""
    if file_data and len(file_data[0]) != len(header[0]):
        raise EnergyFileReaderError(
            f"The header of {filename} names {len(header[0])} columns, "
            f"but the data has {len(file_data[0])}."
        )
