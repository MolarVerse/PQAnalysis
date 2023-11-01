import os

from PQAnalysis.io.base import BaseReader


class EnergyFileReader(BaseReader):
    def __init__(self, filename: str, info_filename: str = None):

        super().__init__(filename)
        self.info_filename = info_filename

        self.withInfoFile = self.__info_file_found__()

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

            self.info_filename = os.path.splitext(self.filename)[0] + ".info"
            try:
                BaseReader(self.info_filename)
            except FileNotFoundError:
                self.info_filename = None
        else:
            try:
                BaseReader(self.info_filename)
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Info File {self.info_filename} not found.")

        if self.info_filename is not None:
            print(f"A Info File \'{self.filename}\' was found.")
            return True
        else:
            return False
