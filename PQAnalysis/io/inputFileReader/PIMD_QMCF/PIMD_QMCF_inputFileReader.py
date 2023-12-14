from ..formats import InputFileFormat
from ..inputFileParser import InputFileParser
from .output_files import _OutputFileMixin


class PIMD_QMCF_InputFileReader(_OutputFileMixin):
    def __init__(self, filename):
        self.format = InputFileFormat.PIMD_QMCF
        self.parser = InputFileParser(filename, self.format)

    def read(self):
        self.dictionary = self.parser.parse()

    def continue_input_file(self, n: int):
        self._parse_old_n()

    def _parse_old_n(self):
        self.old_n = self._get_old_n_from_filename(self.start_file)
        old_n = self._get_old_n_from_filename(self.rpmd_start_file)

    def _get_old_n_from_filename(self, filename: str) -> str:
        filename_splitted = filename.split(".", 1)

        counter = -1
        digits = []
        while isdigit(filename_splitted[counter]):
            digits = [filename_splitted[counter]] + digits
            counter -= 1

        if len(digits) == 0:
            raise ValueError(
                f"Filename {filename} does not contain a number to be continued from. It has to be of the form \"...<number>.<extension>\".")

        return "".join(digits)
