from PQAnalysis.io.base import BaseReader


class InfoFileReader(BaseReader):
    def __init__(self, filename: str):
        """
        Initializes the InfoFileReader with the given filename.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        """
        super().__init__(filename)

    def read(self):
        """
        Reads the info file.

        Returns
        -------
        dict
            The information strings of the info file as a dictionary.
        dict
            The units of the info file as a dictionary.
        """
        info = {}

        with open(self.filename, "r") as file:
            for line in file:
                if line.startswith("--"):
                    continue

                line = line.split()

                if len(line) == 0:
                    continue

                info[line[0]] = line[1:]

        return info
