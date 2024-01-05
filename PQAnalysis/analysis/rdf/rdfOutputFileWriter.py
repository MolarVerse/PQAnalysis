from beartype.typing import Tuple

from . import RadialDistributionFunction
from ...types import Np1DNumberArray
from ...io import BaseWriter


class RDFOutputFileWriter():
    def __init__(self,
                 data_filename: str,
                 log_filename: str | None,
                 data: Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray],
                 rdf: RadialDistributionFunction) -> None:

        self.data_writer = RDFDataWriter(data_filename, data)
        self.log_writer = RDFLogWriter(log_filename, rdf)

    def write(self):
        self.data_writer.write()
        self.log_writer.write()


class RDFDataWriter(BaseWriter):
    def __inti__(self, filename: str,
                 data: Tuple[Np1DNumberArray, Np1DNumberArray,
                             Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]
                 ) -> None:
        self.filename = filename
        super().__init__(filename)

    def write(self):
        pass


class RDFLogWriter(BaseWriter):
    def __inti__(self, filename: str | None, rdf: RadialDistributionFunction) -> None:
        self.filename = filename
        super().__init__(filename)

    def write(self):
        pass
