"""
A module containing the TrajectoryWriter class and its associated methods.
"""

from beartype.typing import List

from PQAnalysis.io.base import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.traj import Trajectory, TrajectoryFormat, MDEngineFormat
from PQAnalysis.core import Cell, Atom
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import runtime_type_checking, runtime_type_checking_setter



class TrajectoryWriter(BaseWriter):

    """
    A class for writing a trajectory to a file.
    Inherits from BaseWriter. See BaseWriter for more information.

    It can write a trajectory to a file in either a PQ format or a QMCFC format.
    """

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
        mode: str | FileWritingMode = 'w'
    ) -> None:
        """
        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If None, the output is printed to stdout.
        engine_format : MDEngineFormat | str, optional
            The format of the md engine for the output file. The default is MDEngineFormat.PQ.
        mode : str, optional
            The mode of the file. Either 'w' for write, 
            'a' for append or 'o' for overwrite. The default is 'w'.
        """

        super().__init__(filename, FileWritingMode(mode))

        self.format = MDEngineFormat(engine_format)

    @runtime_type_checking
    def write(
        self,
        trajectory: Trajectory | AtomicSystem,
        traj_type: TrajectoryFormat | str = TrajectoryFormat.XYZ
    ) -> None:
        """
        Writes the trajectory to the file.

        Parameters
        ----------
        trajectory : Trajectory | AtomicSystem
            The trajectory to write.
        traj_type : TrajectoryFormat | str, optional
            The type of the data to write to the file. Default is TrajectoryFormat.XYZ.
        """

        self.type = TrajectoryFormat(traj_type)

        if isinstance(trajectory, AtomicSystem):
            trajectory = Trajectory([trajectory])

        if self.type == TrajectoryFormat.XYZ:
            self._write_positions(trajectory)
        elif self.type == TrajectoryFormat.VEL:
            self._write_velocities(trajectory)
        elif self.type == TrajectoryFormat.FORCE:
            self._write_forces(trajectory)
        elif self.type == TrajectoryFormat.CHARGE:
            self._write_charges(trajectory)

    def _write_positions(self, trajectory: Trajectory) -> None:
        """
        Writes the positions of the trajectory to the file.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to write.
        """
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.pos, frame.atoms)

        self.close()

    def _write_velocities(self, trajectory: Trajectory) -> None:
        """
        Writes the velocities of the trajectory to the file.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to write.
        """
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.vel, frame.atoms)

        self.close()

    def _write_forces(self, trajectory: Trajectory) -> None:
        """
        Writes the forces of the trajectory to the file.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to write.
        """
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.forces, frame.atoms)

        self.close()

    def _write_charges(self, trajectory: Trajectory) -> None:
        """
        Writes the charges of the trajectory to the file.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to write.
        """
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_scalar(frame.charges, frame.atoms)

        self.close()

    def _write_header(self, n_atoms: int, cell: Cell = Cell()) -> None:
        """
        Writes the header line of the frame to the file.

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the frame.
        cell : Cell
            The cell of the frame. Default is Cell().
        """

        # If the format is QMCFC, an additional atom is added to the count.
        if self.format == MDEngineFormat.QMCFC:
            n_atoms += 1

        if cell != Cell():
            print(
                (
                f"{n_atoms} "
                f"{cell.x} {cell.y} {cell.z} "
                f"{cell.alpha} {cell.beta} {cell.gamma}"
                ),
                file=self.file
            )
        else:
            print(f"{n_atoms}", file=self.file)

    def _write_comment(self, frame: AtomicSystem) -> None:
        """
        Writes the comment line of the frame to the file.

        Parameters
        ----------
        frame : AtomicSystem
            The frame to write the comment line of.
        """

        if self.type == TrajectoryFormat.FORCE:
            sum_forces = sum(frame.forces)
            print(
                (
                f"sum of forces: {sum_forces[0]:e} "
                f"{sum_forces[1]:e} {sum_forces[2]:e}"
                ),
                file=self.file
            )
        else:
            print("", file=self.file)

    def _write_xyz(self, xyz: Np2DNumberArray, atoms: List[Atom]) -> None:
        """
        Writes the xyz of the frame to the file.

        If format is 'qmcfc', an additional X 0.0 0.0 0.0 line is written.

        Parameters
        ----------
        xyz : np.array
            The xyz data of the atoms (either positions, velocities or forces).
        atoms : Elements
            The elements of the frame.
        """

        if self.format == MDEngineFormat.QMCFC and self._type == TrajectoryFormat.XYZ:
            print("X   0.0 0.0 0.0", file=self.file)

        for i, atom in enumerate(atoms):
            if self.type == TrajectoryFormat.VEL:
                print(
                    (
                    f"{atom.name} {xyz[i][0]:16.12e} "
                    f"{xyz[i][1]:16.12e} {xyz[i][2]:16.12e}"
                    ),
                    file=self.file
                )
            else:
                print(
                    (
                    f"{atom.name} {xyz[i][0]:16.10f} "
                    f"{xyz[i][1]:16.10f} {xyz[i][2]:16.10f}"
                    ),
                    file=self.file
                )

    def _write_scalar(
        self,
        scalar: Np1DNumberArray,
        atoms: List[Atom]
    ) -> None:
        """
        Writes the charges of the frame to the file.

        Parameters
        ----------
        scalar : np.array
            scalar data of the atoms (atm only charges).
        atoms : Elements
            The elements of the frame.
        """

        for i, atom in enumerate(atoms):
            print(f"{atom.name} {scalar[i]}", file=self.file)

    @property
    def format(self) -> MDEngineFormat:
        """MDEngineFormat: The format of the md engine for the output file."""
        return self._format

    @format.setter
    @runtime_type_checking_setter
    def format(self, engine_format: MDEngineFormat | str) -> None:
        self._format = MDEngineFormat(engine_format)

    @property
    def type(self) -> TrajectoryFormat:
        """TrajectoryFormat: The type of the data to write to the file."""
        return self._type

    @type.setter
    @runtime_type_checking_setter
    def type(self, traj_type: TrajectoryFormat | str) -> None:
        self._type = TrajectoryFormat(traj_type)
