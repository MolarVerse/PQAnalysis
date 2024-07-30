"""
A module containing the classes for writing related to an
:py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
analysis to a file.
"""

# local imports
from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray
from PQAnalysis.io import BaseWriter
from PQAnalysis.utils import __header__
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.io.formats import FileWritingMode

from .thermal_expansion import ThermalExpansion


class ThermalExpansionDataWriter(BaseWriter):

    """
    Class for writing the data of an
    :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
    analysis to a file.
    """

    @runtime_type_checking
    def __init__(
            self,
            filename: str,
            mode: str | FileWritingMode = "w") -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        mode : str | FileWritingMode, optional
            The writing mode. Default is "w".
            The writing mode can be either a string or a FileWritingMode enum value.
            Possible values are:
            - "w" or FileWritingMode.WRITE: write mode (default, no overwrite)
            - "a" or FileWritingMode.APPEND: append mode
            - "o" or FileWritingMode.OVERWRITE: overwrite mode
        """
        self.filename = filename
        super().__init__(filename, mode=mode)

    @runtime_type_checking
    def write(
        self,
        temperature_points: Np1DNumberArray,
        boxes_avg_data: Np2DNumberArray,
        boxes_std_data: Np2DNumberArray,
        thermal_expansion_data: Np1DNumberArray,

    ):
        """
        Writes the data to the file.

        Parameters
        ----------
        temperature_points : Np1DNumberArray
            the temperature points
        boxes_avg_data : Np2DNumberArray
            the average box data output from the Box._initialize_run() method
        boxes_std_data : Np2DNumberArray
            the standard deviation box data output from the Box._initialize_run() method
        thermal_expansion_data : Np1DNumberArray
            the thermal expansion data output from the Box.run() method
        """

        thermal_expansion_data_mega = thermal_expansion_data * 1e6
        super().open()
        angstrom = '\u212B'.encode('utf-8')
        print(
            f"T / K"
            f"a_avg in {angstrom}      a_std in {angstrom}"
            f"b_avg in {angstrom}      b_std in {angstrom}"
            f"c_avg in {angstrom}      c_std in {angstrom}"
            f"volume in {angstrom}³       volume_std in {angstrom}³"
            f"thermal_expansion_a in (M/K)       thermal_expansion_b in (M/K)"
            f"thermal_expansion_c in (M/K)      volumetric expansion in (M/K)",
            file=self.file
        )
        for i, temperature_point in enumerate(temperature_points):
            print(
                f"{temperature_point}       "
                f"{boxes_avg_data[0][i]}       {boxes_std_data[0][i]}       "
                f"{boxes_avg_data[1][i]}       {boxes_std_data[1][i]}       "
                f"{boxes_avg_data[2][i]}       {boxes_std_data[2][i]}       "
                f"{boxes_avg_data[3][i]}       {boxes_std_data[3][i]}       "
                f"{thermal_expansion_data_mega[0]}       {
                    thermal_expansion_data_mega[1]}       "
                f"{thermal_expansion_data_mega[2]}       {
                    thermal_expansion_data_mega[3]}", file=self.file)
        super().close()


class ThermalExpansionLogWriter(BaseWriter):

    """
    Class for writing the log of an
    :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
    analysis to a file.
    """

    @runtime_type_checking
    def __init__(
            self,
            filename: str,
            mode: str | FileWritingMode = "w") -> None:
        """
        Parameters
        ----------
        filename : str
            the filename to write to
        mode : str | FileWritingMode, optional
            The writing mode. Default is "w".
            The writing mode can be either a string or a FileWritingMode enum value.
            Possible values are:
            - "w" or FileWritingMode.WRITE: write mode (default, no overwrite)
            - "a" or FileWritingMode.APPEND: append mode
            - "o" or FileWritingMode.OVERWRITE: overwrite mode
        """
        self.filename = filename
        super().__init__(filename, mode=mode)

    @runtime_type_checking
    def write_before_run(self, thermal_expansion: ThermalExpansion) -> None:
        """
        Writes the log before the
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
        run() method is called.

        This includes the general header of PQAnalysis
        and the most important setup parameters of the
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
        analysis.

        Parameters
        ----------
        thermal_expansion : ThermalExpansion
            the thermal expansion coefficient object
        """
        super().open()

        if self.filename is not None:
            print(__header__, file=self.file)
            print(file=self.file)

        print("Thermal Expansion Analysis", file=self.file)

        print(f"    Number of temperature points: {
              len(thermal_expansion.temperature_points)}", file=self.file)
        print(f"    Temperature points: {
            thermal_expansion.temperature_points}", file=self.file)
        print(f"    Temperature step size: {
            thermal_expansion.temperature_step_size}", file=self.file)
        print(f"    Number of box files: {
            len(thermal_expansion.boxes)}", file=self.file)

        print(file=self.file)

        super().close()

    @runtime_type_checking
    def write_after_run(self, thermal_expansion: ThermalExpansion) -> None:
        """
        Writes the log after the
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
        run() method is called.

        This includes the elapsed time of the
        :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`
        run() method.

        Parameters
        ----------
        thermal_expansion : ThermalExpansion
            the linear thermal expansion coefficient object
        """
        super().open()
        angstrom = '\u212B'.encode('utf-8')
        print(
            f"    Average a: {thermal_expansion.boxes_avg[0]}{
                angstrom}+/- {thermal_expansion.boxes_std[0]}{angstrom}",
            file=self.file
        )
        print(
            f"    Average b: {thermal_expansion.boxes_avg[1]}{
                angstrom}+/- {thermal_expansion.boxes_std[1]}{angstrom}",
            file=self.file
        )
        print(
            f"    Average c: {thermal_expansion.boxes_avg[2]}{
                angstrom}+/- {thermal_expansion.boxes_std[2]}{angstrom}",
            file=self.file
        )
        print(
            f"    Average volume: {thermal_expansion.boxes_avg[3]}{
                angstrom}³ +/- {thermal_expansion.boxes_std[3]}{angstrom}³",
            file=self.file
        )
        print(
            f"    Middle points: {thermal_expansion.middle_points} at {
                thermal_expansion.temperature_points[
                    len(thermal_expansion.temperature_points) // 2
                ]} K",
            file=self.file
        )
        print(
            f"    Linear thermal expansion a: {
                thermal_expansion.thermal_expansions[0]} 1/K",
            file=self.file
        )
        print(
            f"    Linear thermal expansion b: {
                thermal_expansion.thermal_expansions[1]} 1/K",
            file=self.file
        )
        print(
            f"    Linear thermal expansion c: {
                thermal_expansion.thermal_expansions[2]} 1/K",
            file=self.file
        )
        print(
            f"    Volumetric expansion: {
                thermal_expansion.thermal_expansions[3]} 1/K",
            file=self.file
        )

        print(file=self.file)

        print(
            "    Linear thermal expansion a in M/K: ",
            file=self.file
        )
        print(
            f"    {thermal_expansion.thermal_expansions[0] * 1e6} M/K",
            file=self.file
        )
        print(
            "    Linear thermal expansion b in M/K: ", file=self.file
        )
        print(
            f"    {thermal_expansion.thermal_expansions[1] * 1e6} M/K",
            file=self.file
        )
        print(
            "    Linear thermal expansion c in M/K ", file=self.file
        )
        print(
            f"    {thermal_expansion.thermal_expansions[2] * 1e6} M/K",
            file=self.file
        )

        print(
            "    Volumetric thermal expansion in M/K: ",
            file=self.file
        )
        print(
            f"    {thermal_expansion.thermal_expansions[3] * 1e6} M/K",
            file=self.file
        )

        self.file.write(f"Elapsed time: {
            thermal_expansion.elapsed_time} s\n")

        super().close()
