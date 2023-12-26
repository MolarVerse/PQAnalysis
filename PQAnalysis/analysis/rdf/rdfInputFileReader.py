from ...io import InputFileParser
from ...io.inputFileReader.inputFileParser import InputDictionary
from ...io.inputFileReader.exceptions import InputFileError
from ...io.inputFileReader.PQAnalysis.parse import parse_traj_files_key


class RDFInputFileReader:

    required_keys = ["traj_files", "reference_selection", "target_selection"]

    def __init__(self, filename: str) -> None:
        self.filename = filename
        self.parser = InputFileParser(filename)

    def read(self):
        self.dictionary = self.parser.parse()

        if not all([key in self.dictionary.keys() for key in self.required_keys]):
            raise InputFileError(
                f"Not all required keys were set in the input file! The required keys are: {self.required_keys}.")

        self.traj_files = parse_traj_files_key(self.dictionary)
