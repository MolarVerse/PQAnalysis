from beartype.typing import Tuple

from . import RadialDistributionFunction
from ...types import Np1DNumberArray
from ...io import BaseWriter
from ...utils import header


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
        self.data = data
        super().__init__(filename)

    def write(self):
        super().open()

        for i in range(len(self.data[0])):
            print(f"{self.data[0][i]} {self.data[1][i]} {self.data[2][i]} {self.data[3][i]} {self.data[4][i]}",
                  file=self.file)

        super().close()


class RDFLogWriter(BaseWriter):
    def __inti__(self, filename: str | None, rdf: RadialDistributionFunction) -> None:
        self.filename = filename
        super().__init__(filename)

    def write_before_run(self):
        super().open()

        print(header, file=self.file)
        print(file=self.file)

        print("RDF calculation:", file=self.file)
        print(file=self.file)

        angstrom = u'\u212B'.encode('utf-8')

        # fmt: off
        print(f"    Number of bins: {self.rdf.n_bins}", file=self.file)
        print(f"    Bin width:      {self.rdf.delta_r} {angstrom}", file=self.file)
        print(f"    Minimum radius: {self.rdf.r_min} {angstrom}", file=self.file)
        print(f"    Maximum radius: {self.rdf.r_max} {angstrom}", file=self.file)
        print(file=self.file)
        # fmt: on

        print(f"    Number of frames: {self.rdf.n_frames}", file=self.file)
        print(f"    Number of atoms:  {self.rdf.n_atoms}", file=self.file)
        print(file=self.file)

        # fmt: off
        print(f"    Reference selection: {self.rdf.reference_selection}", file=self.file)
        print(f"    total number of atoms in reference selection: {len(self.rdf.reference_indices)}", file=self.file)
        print(f"    Target selection:    {self.rdf.target_selection}", file=self.file)
        print(f"    total number of atoms in target selection:    {len(self.rdf.target_indices)}", file=self.file)
        # fmt: on

        super().close()

    def write_after_run(self):
        super().open()

        print(f"    Elapsed time: {self.rdf.elapsed_time} ms", file=self.file)

        super().close()
