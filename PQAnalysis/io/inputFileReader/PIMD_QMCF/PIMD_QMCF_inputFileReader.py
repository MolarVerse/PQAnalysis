from ..formats import InputFileFormat
from ..inputFileParser import InputFileParser
from .output_files import _OutputFileMixin


class PIMD_QMCF_InputFileReader(_OutputFileMixin):
    def __init__(self, filename):
        self.format = InputFileFormat.PIMD_QMCF
        self.parser = InputFileParser(filename, self.format)

    def read(self):
        self.dictionary = self.parser.parse()
