"""
A module containing the PIMD_QMCF_InputFileReader class.

Classes
-------
PIMD_QMCF_InputFileReader
    Reads a PIMD_QMCF input file and parses it.
    
Functions
---------
_increase_digit_string
    increases a string containing only digits by one
_get_digit_string_from_filename
    extracts a string containing only digits from a filename
"""

import re

from ....types import PositiveInt
from ..formats import InputFileFormat
from ..inputFileParser import InputFileParser
from .output_files import _OutputFileMixin


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
    ValueError
        if digit_string contains non-digit characters
    """

    if not all([char.isdigit() for char in digit_string]):
        raise ValueError(
            f"digit_string {digit_string} contains non-digit characters.")

    string_without_leading_zeros = digit_string.lstrip("0")

    if string_without_leading_zeros == "":
        string_without_leading_zeros = "0"

    string_without_leading_zeros = str(int(string_without_leading_zeros) + 1)

    return "0" * (len(digit_string) - len(string_without_leading_zeros)) + string_without_leading_zeros


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
    "filename.extension" -> ValueError

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
    ValueError
        if filename does not contain a number to be parsed
    """

    regex = re.search(r"\d+.", filename)

    if regex is None:
        raise ValueError(
            f"Filename {filename} does not contain a number to be continued from. It has to be of the form \"...<number>.<extension>\".")

    return regex.group(0)[:-1]


class PIMD_QMCF_InputFileReader(_OutputFileMixin):
    """
    Reads a PIMD_QMCF input file and parses it.

    Parameters
    ----------
    _OutputFileMixin : class
        mixin class containing all output and start file keys
    """

    def __init__(self, filename: str):
        """
        Initialize the PIMD_QMCF_InputFileReader class.

        self.format is set to InputFileFormat.PIMD_QMCF.
        A parser is created using the InputFileParser class.

        Parameters
        ----------
        filename : str
            filename of the input file
        """

        self.format = InputFileFormat.PIMD_QMCF
        self.filename = filename
        self.parser = InputFileParser(self.filename, self.format)

    def read(self):
        """
        Reads the input file and parses it.

        The dictionary is set to the dictionary returned by the parser.
        The raw_input_file is set to the raw_input_file returned by the parser.
        """

        self.dictionary = self.parser.parse()
        self.raw_input_file = self.parser.raw_input_file

        if not self.is_start_file_defined:
            raise ValueError(
                f"No start file defined in input file {self.filename}.")

    def continue_input_file(self, n: PositiveInt):
        """
        Creates n new input files by increasing the number in the filename by one.
        All other numbers in the start- and output-files within the input file are increased by one as well.

        Parameters
        ----------
        n : int
            number of new input files to be created

        Raises
        ------
        ValueError
            if the n parsed from the output files defined in the input file does not match the n parsed from the input file name
        ValueError
            if the n parsed from the start file does not match the n parsed from the output files
        """
        self.input_file_n = _get_digit_string_from_filename(self.filename)
        self.start_n = self._parse_start_n()
        self.actual_n = self._parse_actual_n()

        # check if n from input file name matches n from output files
        if int(self.actual_n) != int(self.input_file_n):
            raise ValueError(
                f"Actual n ({self.actual_n}) and input file n ({self.input_file_n}) do not match.")

        # check if n from start file matches n from output files
        if int(self.start_n) != int(self.actual_n) - 1:
            raise ValueError(
                f"Old n ({self.start_n}) has to be one less than actual n ({self.actual_n}).")

        old_input_file_n = self.input_file_n
        old_start_n = self.start_n
        old_actual_n = self.actual_n

        for _ in range(n):
            new_input_file_n = _increase_digit_string(old_input_file_n)
            new_start_n = _increase_digit_string(old_start_n)
            new_actual_n = _increase_digit_string(old_actual_n)

            new_raw_input_file = self.raw_input_file

            # replace digit strings in start files
            for key in self.output_file_keys:
                if key in self.dictionary.keys():

                    new_filename = self.dictionary[key][0].replace(
                        self.actual_n, new_actual_n)

                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename)

            # replace digit strings in output files
            for key in self.start_file_keys:
                if key in self.dictionary.keys():

                    new_filename = self.dictionary[key][0].replace(
                        self.start_n, new_start_n)

                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename)

            # create new input file and write new_raw_input_file to it
            new_filename = self.filename.replace(
                self.input_file_n, new_input_file_n)
            file = open(new_filename, "w")
            file.write(new_raw_input_file)
            file.close()

            old_input_file_n = new_input_file_n
            old_start_n = new_start_n
            old_actual_n = new_actual_n

    def _parse_start_n(self) -> str:
        """
        Parses the n from the start file.

        If the rpmd_start_file is defined, the n from the start_file and the rpmd_start_file are compared.
        If they do not match, a ValueError is raised.

        Returns
        -------
        str
            n parsed from the start file as a digit string

        Raises
        ------
        ValueError
            if the n from the start file and the rpmd_start_file do not match
        """
        n = _get_digit_string_from_filename(self.start_file)

        if self.is_rpmd_start_file_defined:
            # add "." to match also files without extension
            _n = _get_digit_string_from_filename(self.rpmd_start_file + ".")

            if _n != n:
                raise ValueError(
                    f"N from start_file ({n}) and rpmd_start_file ({_n}) do not match.")

        return n

    def _parse_actual_n(self) -> str:
        """
        Parses the n from the output files.

        Returns
        -------
        str
            n parsed from the output files as a digit string

        Raises
        ------
        ValueError
            if no output file is defined
        ValueError
            if the n parsed from the output files is not consistent
        """
        n = None

        for key in self.output_file_keys:
            if key in self.dictionary.keys():
                # add "." to match also files without extension
                _n = _get_digit_string_from_filename(
                    self.dictionary[key][0] + ".")

                if _n != n and n is not None:
                    print(_n, n, key)
                    raise ValueError(
                        f"Actual n in output files is not consistent.")

                n = _n

        if n is None:
            raise ValueError(
                f"No output file found to determine actual n.")

        return n
