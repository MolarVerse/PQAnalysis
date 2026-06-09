"""
A module containing the PQ_InputFileReader class.
"""

import logging
import re

from PQAnalysis.types import PositiveInt
from PQAnalysis.io.input_file_reader.formats import InputFileFormat
from PQAnalysis.io.input_file_reader.input_file_parser import InputFileParser
from PQAnalysis.exceptions import PQNotImplementedError, PQValueError
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from .output_files import _OutputFileMixin



class PQInputFileReader(_OutputFileMixin):

    """
    Reads a PQ input file and parses it.

    Parameters
    ----------
    _OutputFileMixin : class
        mixin class containing all output and start file keys
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    def __init__(
        self,
        filename: str,
        input_format: InputFileFormat | str = InputFileFormat.PQ
    ):
        """
        Initialize the PQ_InputFileReader class.

        self.format is set to a QMCF input format.
        A parser is created using the InputFileParser class.

        Parameters
        ----------
        filename : str
            filename of the input file
        input_format : InputFileFormat | str, optional
            the format of the input file, by default InputFileFormat.PQ
        """

        ########################
        # dummy initialization #
        ########################

        self.dictionary = None
        self.raw_input_file = None
        self.input_file_n = None
        self.start_n = None
        self.actual_n = None

        self.format = InputFileFormat(input_format)

        if not InputFileFormat.is_qmcf_type(self.format):
            self.logger.error(
                (
                    f"Format {self.format} not implemented "
                    "for PQ input file reading."
                ),
                exception=PQNotImplementedError
            )

        self.filename = filename
        self.parser = InputFileParser(self.filename, self.format)

    def read(self):
        """
        Reads the input file and parses it.

        The dictionary is set to the dictionary returned by the parser.
        The raw_input_file is set to the raw_input_file returned by the parser.

        Raises
        ------
        PQValueError
            if no start file is defined in the input file
        """

        self.dictionary = self.parser.parse()
        self.raw_input_file = self.parser.raw_input_file

        if not self.is_start_file_defined:
            self.logger.error(
                f"No start file defined in input file {self.filename}.",
                exception=PQValueError
            )

    def continue_input_file(self, n: PositiveInt):
        """
        Creates n new input files by increasing the number in the 
        filename by one. All other numbers in the start- and
        output-files within the input file are increased by one as well.

        Parameters
        ----------
        n : int
            number of new input files to be created

        Raises
        ------
        PQValueError
            if the n parsed from the output files defined in the input file does
            not match the n parsed from the input file name
        PQValueError
            if the n parsed from the start file does not match the n parsed from the output files
        """
        self.input_file_n = _get_digit_string_from_filename(self.filename)
        start_file_ns = self._parse_start_file_ns()
        self.start_n = next(
            (
                start_n for start_n in start_file_ns.values()
                if start_n is not None
            ),
            None
        )
        self.actual_n = self._parse_actual_n()

        # check if n from input file name matches n from output files
        if int(self.actual_n) != int(self.input_file_n):
            self.logger.error(
                f"Actual n ({self.actual_n}) and input file n ({self.input_file_n}) do not match.",
                exception=PQValueError
            )

        # check if n from start file matches n from output files
        if self.start_n is not None and int(self.start_n) != int(self.actual_n) - 1:
            self.logger.error(
                f"Old n ({self.start_n}) has to be one less than actual n ({self.actual_n}).",
                exception=PQValueError
            )

        old_input_file_n = self.input_file_n
        old_start_n = self.start_n
        old_actual_n = self.actual_n

        for _ in range(n):
            new_input_file_n = _increase_digit_string(old_input_file_n)
            new_start_n = (
                _increase_digit_string(old_start_n)
                if old_start_n is not None
                else old_actual_n
            )
            new_actual_n = _increase_digit_string(old_actual_n)

            new_raw_input_file = self.raw_input_file

            # replace digit strings in start files
            for key in self.output_file_keys:
                if key in self.dictionary.keys():

                    new_filename = self.dictionary[key][0].replace(
                        self.actual_n, new_actual_n
                    )

                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename
                    )

            # replace digit strings in output files
            for key in self.start_file_keys:
                if key in self.dictionary.keys():

                    if start_file_ns[key] is None:
                        new_filename = self._continued_start_file(key, new_start_n)
                    else:
                        new_filename = self.dictionary[key][0].replace(
                            start_file_ns[key], new_start_n
                        )

                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename
                    )

            # create new input file and write new_raw_input_file to it
            new_filename = self.filename.replace(
                self.input_file_n, new_input_file_n
            )

            with open(new_filename, "w", encoding="utf-8") as file:
                file.write(new_raw_input_file)
                file.close()

            old_input_file_n = new_input_file_n
            old_start_n = new_start_n
            old_actual_n = new_actual_n

    def _parse_start_file_ns(self) -> dict[str, str | None]:
        """
        Parses the n from the start files.

        If multiple start files define a number, the numbers are compared. If
        they do not match, a PQValueError is raised. Start files without a
        parseable number are allowed and are continued from restart output.

        Returns
        -------
        dict[str, str | None]
            n parsed from each start file as a digit string or None if the
            filename has no parseable digit string

        Raises
        ------
        PQValueError
            if the n from multiple start files do not match
        """
        start_file_ns = {}

        for key in self.start_file_keys:
            if key in self.dictionary.keys():
                start_file_ns[key] = _get_optional_digit_string_from_filename(
                    self.dictionary[key][0] + "."
                )

        n = next(
            (start_n for start_n in start_file_ns.values() if start_n is not None),
            None
        )

        for key, start_n in start_file_ns.items():
            if start_n is not None and start_n != n:
                self.logger.error(
                    (
                        f"N from start_file ({n}) and "
                        f"{key} ({start_n}) do not match."
                    ),
                    exception=PQValueError
                )

        return start_file_ns

    def _parse_start_n(self) -> str:
        """
        Parses the n from the start file.

        If the rpmd_start_file is defined, the n from the start_file and 
        the rpmd_start_file are compared. If they do not match,
        a PQValueError is raised.

        Returns
        -------
        str
            n parsed from the start file as a digit string

        Raises
        ------
        PQValueError
            if the n from the start file and the rpmd_start_file do not match
        """
        n = _get_digit_string_from_filename(self.start_file)

        if self.is_rpmd_start_file_defined:
            # add "." to match also files without extension
            if (
                _n :=
                _get_digit_string_from_filename(self.rpmd_start_file + '.')
            ) != n:
                self.logger.error(
                    f"N from start_file ({n}) and rpmd_start_file ({_n}) do not match.",
                    exception=PQValueError
                )

        return n

    def _continued_start_file(self, key: str, n: str) -> str:
        """
        Gets the next start file name for an unnumbered start file.
        """
        restart_key = "restart_file"

        if key == "rpmd_start_file":
            restart_key = "rpmd_restart_file"

        if restart_key not in self.dictionary.keys() and key == "start_file":
            restart_key = "rpmd_restart_file"

        if restart_key not in self.dictionary.keys():
            if "file_prefix" not in self.dictionary.keys():
                return self.dictionary[key][0]

            file_prefix = self.dictionary["file_prefix"][0].replace(
                self.actual_n,
                n
            )

            if key == "rpmd_start_file":
                return f"{file_prefix}.rpmd.rst"

            return f"{file_prefix}.rst"

        return self.dictionary[restart_key][0].replace(self.actual_n, n)

    def _parse_actual_n(self) -> str:
        """
        Parses the n from the output files.

        Returns
        -------
        str
            n parsed from the output files as a digit string

        Raises
        ------
        PQValueError
            if no output file is defined
        PQValueError
            if the n parsed from the output files is not consistent
        """
        n = None

        for key in self.output_file_keys:
            if key in self.dictionary.keys():
                # add "." to match also files without extension
                _n = _get_digit_string_from_filename(
                    self.dictionary[key][0] + "."
                )

                if _n != n and n is not None:
                    self.logger.error(
                        "Actual n in output files is not consistent.",
                        exception=PQValueError
                    )

                n = _n

        if n is None:
            self.logger.error(
                "No output file found to determine actual n.",
                exception=PQValueError
            )

        return n



def _increase_digit_string(digit_string: str) -> str:
    """
    increases a string containing only digits by one

    If the string contains leading zeros, they are preserved.
    The total length of the string is preserved, if there are leading zeros.

    For example:
    "0001" -> "0002"
    "001" -> "002"
    "099" -> "100"
    "999" -> "1000"

    Parameters
    ----------
    digit_string : str
        string containing only digits to be increased

    Returns
    -------
    str
        increased string by one

    Raises
    ------
    PQValueError
        if digit_string contains non-digit characters
    """

    if not all(char.isdigit() for char in digit_string):
        PQInputFileReader.logger.error(
            f"digit_string {digit_string} contains non-digit characters.",
            exception=PQValueError
        )

    if (without_leading_zeros := digit_string.lstrip('0')) == '':
        without_leading_zeros = "0"

    without_leading_zeros = str(int(without_leading_zeros) + 1)

    number_leading_zeros = len(digit_string) - len(without_leading_zeros)
    leading_zeros = "0" * number_leading_zeros

    return leading_zeros + without_leading_zeros



def _get_digit_string_from_filename(filename: str) -> str:
    """
    extracts a string containing only digits from a filename

    The string has to be of the form "...<number>.<extension>".
    The string is extracted from the number part of the filename.

    For example:
    "filename_0001.extension" -> "0001"
    "filename_001.extension" -> "001"
    "filename_099.extension" -> "099"
    "filename_100.extension" -> "100"
    "filename.extension" -> PQValueError

    Parameters
    ----------
    filename : str
        filename to be parsed

    Returns
    -------
    str
        string containing only digits

    Raises
    ------
    PQValueError
        if filename does not contain a number to be parsed
    """

    if (digit_string := _get_optional_digit_string_from_filename(filename)) is None:
        PQInputFileReader.logger.error(
            (
                f"Filename {filename} does not contain a number to be "
                "continued from. It has to be of the form \"...<number>.<extension>\"."
            ),
            exception=PQValueError
        )

    return digit_string



def _get_optional_digit_string_from_filename(filename: str) -> str | None:
    """
    Extracts a digit string from a filename if one can be parsed.
    """
    if (regex := re.search(r"\d+\.", filename)) is None:
        return None

    return regex.group(0)[:-1]
