from ..formats import InputFileFormat
from ..inputFileParser import InputFileParser
from .output_files import _OutputFileMixin


def _continue_digit_string(digit_string: str) -> str:
    n_length = len(digit_string)
    leading_zeros = 0

    for char in digit_string:
        if char != "0":
            break
        else:
            leading_zeros += 1

    digit_string = str(int(digit_string[leading_zeros:]) + 1)

    delta_n = n_length - len(digit_string)
    for _ in range(delta_n):
        digit_string = "0" + digit_string

    return digit_string


def _get_n_from_filename(filename: str) -> str:
    filename_splitted = filename.split(".", 1)[0]

    counter = -1
    digits = []
    while filename_splitted[counter].isdigit():
        digits = [filename_splitted[counter]] + digits
        counter -= 1

    if len(digits) == 0:
        raise ValueError(
            f"Filename {filename} does not contain a number to be continued from. It has to be of the form \"...<number>.<extension>\".")

    n = "".join(digits)

    return n


class PIMD_QMCF_InputFileReader(_OutputFileMixin):
    def __init__(self, filename):
        self.format = InputFileFormat.PIMD_QMCF
        self.filename = filename
        self.parser = InputFileParser(self.filename, self.format)

    def read(self):
        self.dictionary = self.parser.parse()
        self.raw_input_file = self.parser.raw_input_file

    def continue_input_file(self, n: int):
        self.input_file_n = _get_n_from_filename(self.filename)
        self.old_n = self._parse_old_n()
        self.actual_n = self._parse_actual_n()

        if int(self.actual_n) != int(self.input_file_n):
            raise ValueError(
                f"Actual n ({self.actual_n}) and input file n ({self.input_file_n}) do not match.")

        if int(self.old_n) != int(self.actual_n) - 1:
            raise ValueError(
                f"Old n ({self.old_n}) has to be one less than actual n ({self.actual_n}).")

        if int(self.input_file_n) >= n:
            raise ValueError(
                f"Input file n ({self.input_file_n}) has to be less than n ({n}).")

        for _ in range(int(self.actual_n), n + 1):
            new_input_file_n = _continue_digit_string(self.input_file_n)
            new_old_n = _continue_digit_string(self.old_n)
            new_actual_n = _continue_digit_string(self.actual_n)

            new_raw_input_file = self.raw_input_file

            new_filename = self.filename.replace(
                self.input_file_n, new_input_file_n)
            file = open(new_filename, "w")

            for key in self.start_file_keys:
                if key in self.dictionary.keys():
                    new_filename = self.dictionary[key][0].replace(
                        self.old_n, new_old_n)
                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename)

            for key in self.output_file_keys:
                if key in self.dictionary.keys():
                    new_filename = self.dictionary[key][0].replace(
                        self.actual_n, new_actual_n)
                    new_raw_input_file = new_raw_input_file.replace(
                        self.dictionary[key][0], new_filename)

            file.write(new_raw_input_file)
            file.close()

            self.input_file_n = new_input_file_n
            self.old_n = new_old_n
            self.actual_n = new_actual_n

    def _parse_old_n(self) -> str:
        old_n = _get_n_from_filename(self.start_file)

        if self.is_rpmd_start_file_defined:
            _old_n = _get_n_from_filename(self.rpmd_start_file)

            if _old_n != old_n:
                raise ValueError(
                    f"Old n from start_file ({self.old_n}) and rpmd_start_file ({old_n}) do not match.")

        return old_n

    def _parse_actual_n(self) -> str:
        actual_n = None

        for key in self.output_file_keys:
            if key in self.dictionary.keys():
                _actual_n = _get_n_from_filename(self.dictionary[key][0])

                if _actual_n != actual_n and actual_n is not None:
                    raise ValueError(
                        f"Actual n in output files is not consistent.")

                actual_n = _actual_n

        if actual_n is None:
            raise ValueError(
                f"No output file found to determine actual n.")

        return actual_n
