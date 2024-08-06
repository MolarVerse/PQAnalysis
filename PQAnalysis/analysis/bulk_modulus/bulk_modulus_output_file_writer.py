"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` analysis to a file.
"""

# 3rd party imports
from beartype.typing import Tuple

# local imports
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking

from .bulk_modulus import BulkModulus


class BulkModulusDataWriter(BaseWriter):

    """
    Class for writing the data of an 
    :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`
    analysis to a file.
    """

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        """
        self.filename = filename
        super().__init__(filename)

    @runtime_type_checking
    def write(
        self,
        data: Tuple[Np1DNumberArray,
                    Np1DNumberArray,
                    Np1DNumberArray,
                    ]
    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        data : Tuple[Np1DNumberArray, Np1DNumberArray,
            Np1DNumberArray]
            the data output from the BulkModulus.run() method
        """
        super().open()

        for i in range(len(data[0])):
            print(
                (
                    f"{data[0][i]} {data[1][i]} {data[2][i]} "
                ),
                file=self.file
            )

        super().close()


class BulkModulusLogWriter(BaseWriter):

    """
    Class for writing the log (setup parameters) of an 
    :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` analysis
    to a file.
    """

    @runtime_type_checking
    def __init__(self, filename: str | None) -> None:
        """
        Parameters
        ----------
        filename : str | None
            the filename to write to if None, the output is printed to stdout
        """
        self.filename = filename
        super().__init__(filename)

    @runtime_type_checking
    def write_before_run(self, bulk_modulus: BulkModulus):
        """
        Writes the log before the 
        :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` analysis.

        Parameters
        ----------
        bulk_modulus : BulkModulus
            the BulkModulus analysis object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("Bulk Modulus calculation:", file=self.file)
        print(file=self.file)

        print("Setup parameters:", file=self.file)
        print(file=self.file)

        print(f"    Mode: {bulk_modulus.mode}", file=self.file)

        angstrom = '\u212B'.encode('utf-8')

        print(f"    Equilibrium volume: {bulk_modulus.volume_equilibrium} {angstrom}^3",
              file=self.file)
        for i in range(len(bulk_modulus.volumes_perturbation)):
            print(f"    Perturbation volume {i}: {bulk_modulus.volumes_perturbation[i]} {angstrom}^3",
                  file=self.file)




        # fmt: off


        super().close()

    @runtime_type_checking
    def write_after_run(self, bulk_modulus: BulkModulus, units: str = "bar"):
        """
        Writes the log after the 
        :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`
        run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`
        run() method.

        Parameters
        ----------
        bulk_modulus : BulkModulus
            the BulkModulus analysis object
        """
        super().open()

        print(file=self.file)

        for i in range(len(bulk_modulus.pressures_perturbation_avg)):
            print(f"    Perturbation pressure {i}: {bulk_modulus.pressures_perturbation_avg[i]} +/- {bulk_modulus.pressures_perturbation_std[i]} {units}",
                  file=self.file)
        # bulk modulus results
        print("Bulk Modulus results:", file=self.file)
        print(file=self.file)

        print(f"    Bulk modulus: {bulk_modulus.bulk_modulus} +/- {bulk_modulus.bulk_modulus_err} {units}",
              file=self.file)

        if bulk_modulus.mode == "bulk_modulusEOS":
            print(f"    Bulk modulus prime: {bulk_modulus.bulk_modulus_prime} +/- {bulk_modulus.bulk_modulus_prime_err} {units}",
                  file=self.file)

        print(file=self.file)




        print(f"    Elapsed time: {
              bulk_modulus.elapsed_time} s", file=self.file)

        super().close()
