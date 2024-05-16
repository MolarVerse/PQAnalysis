"""
A module containing the MSD class. 
"""

import logging

# 3rd party imports
import numpy as np

# 3rd party imports
from beartype.typing import Tuple
from tqdm.auto import tqdm


# local absolute imports
from PQAnalysis.config import with_progress_bar
from PQAnalysis.types import Np2DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cells
from PQAnalysis.traj import Trajectory, check_trajectory_pbc, check_trajectory_vacuum 
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import MSDError


class MSD:
    """
    A class for calculating the Mean Squared Displacement (MSD) of a trajectory.
    """
    _frame_start_default = 0
    _use_full_atom_default = False
    _no_intra_molecular_default = False
    _window_size_default = 2500
    _gap_default = 10

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self,
                 traj: Trajectory | TrajectoryReader,
                 target_species: SelectionCompatible,
                 use_full_atom_info: bool = False,
                 no_intra_molecular: bool = False,
                 window_size: PositiveInt | None = 500,
                 gap: PositiveInt | None = 1,
                 frame_start: PositiveInt | None = None,
                 frame_stop: PositiveInt | None = None,
                 n_frames: PositiveInt | None = None,
                 ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory to analyze. If a TrajectoryReader is provided,
            the trajectory frame by frame via a frame_generator
        target_species : SelectionCompatible
            The target species of the MSD analysis.
        use_full_atom_info : bool, optional
            Whether to use the full atom information of the trajectory
            or not, by default None (False).
        no_intra_molecular : bool, optional
            Whether to exclude intra-molecular distances or not, by default None (False).
        window_size : PositiveInt, optional
            The window size of the running average, by default 500.
        gap : PositiveInt, optional
            The gap of the running average, by default 1.
        frame_start : PositiveInt, optional
            The starting frame of the trajectory, by default None.
        frame_stop : PositiveInt, optional
            The stopping frame of the trajectory, by default None.
        n_frames : PositiveInt, optional
            The number of frames to use, by default None.


        Raises
        ------
        MSDError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.
        MSDError
            If the trajectory is empty.

        Notes
        -----        
        Furthermore, 

        See Also
        --------
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        :py:class:`~PQAnalysis.topology.selection.Selection`
        :py:class:`~PQAnalysis.io.trajectoryReader.TrajectoryReader`
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        """

        ##############
        # dummy init #
        ##############

        self._origin_positions = np.array([])

        self._window = Trajectory()

        self.msd = np.array([])

        #####################################################
        # Initialize parameters with default values if None #
        #####################################################

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        if no_intra_molecular is None:
            self.no_intra_molecular = self._no_intra_molecular_default
        else:
            self.no_intra_molecular = no_intra_molecular
        
        if window_size is None:
            self.window_size = self._window_size_default
        else:
            self.window_size = window_size
        
        if gap is None:
            self.gap = self._gap_default
        else:
            self.gap = gap

        if n_frames is None and frame_stop is None:
            raise MSDError(
                "Either n_frames or frame_stop has to be defined.")
        
        elif n_frames is not None and frame_stop is not None:
            raise MSDError(
                "Either n_frames or frame_stop has to be defined, not both.")
        
        elif n_frames is not None:
            self.n_frames = n_frames
            self.frame_start = self._frame_start_default
            self.frame_stop = self.frame_start + self.n_frames - 1

        elif frame_stop is not None:
            self.frame_stop = frame_stop
            self.n_frames = self.frame_stop - self.frame_start
    
        ################################
        # Initialize Selection objects #
        ################################

        self.target_species = target_species

        self.target_selection = Selection(target_species)

        ############################################
        # Initialize Trajectory iterator/generator #
        ############################################

        self.cells = traj.cells

        if isinstance(traj, TrajectoryReader):

            self.frame_generator = traj.frame_generator()
        elif len(traj) > 0:
            # use trajectory object as iterator
            self.frame_generator = iter(traj)
        else:
            self.logger.error(
                "Trajectory cannot be of length 0.",
                exception=MSDError
            )
        self.first_frame = next(self.frame_generator)
        self.topology = traj.topology
        self.target_indices = self.target_selection.select(
            self.topology, self.use_full_atom_info)
        
        self._number_windows = self.n_frames - self.window_size + 1
        """
        Sets up the bins of the MSD analysis.

        This method is called by the __init__ method of the 
    
        Parameters
        ----------
        

        Raises
        ------
        MSDError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.

        """


        self._check_trajectory_conditions()


    def _check_trajectory_conditions(self):
        """
        Checks if the trajectory is fully periodic or fully in vacuum.

        Raises
        ------
        MSDError
            If the trajectory is not fully periodic or fully in vacuum.
            Meaning that some frames are in vacuum and others are periodic.
        """

        if not check_trajectory_pbc(self.cells) and not check_trajectory_vacuum(self.cells):
            self.logger.error(
                (
                    "The provided trajectory is not fully periodic or "
                    "in vacuum, meaning that some frames are in vacuum "
                    "and others are periodic. This is not supported by "
                    "the MSD analysis."
                ),
                exception=MSDError
            )

    @timeit_in_class
    def run(self) -> Tuple[
        Np2DNumberArray,
    ]:
        """
        Runs the MSD analysis.



        Returns
        -------

        """

        self._initialize_run()
        self._run()
        return self._finalize_run()
   
            
    
#    def _unwrap_positions(self,method:str="hlat") -> Trajectory:
#         """
#         Unwraps the positions of the target species in the window trajectory

#         """
#         if method.lower() == "hlat":
#             return self._window.unwrap_positions_hlat()
#         else:
#             self.logger.error(
#                 "Method not supported for unwrapping positions.",
#                 exception=MSDError
# #             )
    def _calculate_msd(self, window_start: int = 0, window_end: int | None = None) -> None:
            """
            Generates a window of frames for the MSD analysis.
            """
            for i, frame in enumerate(self.frame_generator, start=window_start, end=window_end):

                if i == 0:
                    self._window = frame.pos[self.target_indices]
                    continue

                if i % self.gap != 0:
                    continue

                if i % self.window_size == 0:
                    self._calculate_msd()
                    self._window = frame.pos[self.target_indices]
                    continue
                unwrap_pos = frame.pos - frame.cell.image(frame.pos[self.target_indices] - self._window[i - 1].pos)
                self._window = self._window.append(frame.replace(pos=unwrap_pos))   
                distances = distance(self._origin_positions, self._window[i].pos)
                self.msd = np.append(self.msd, np.mean(distances**2))

      


    def _initialize_run(self):
        """
        Initializes the MSD analysis for running.

        """
        self._origin_positions = self.first_frame.pos[self.target_indices]


        self._calculate_msd(0, self.window_size)



     
    def _run(self):
        """
        Calculates the bins of the MSD analysis.

   
        """
       
        for i in range(1, self._number_windows):
            self._calculate_msd(i, i + self.window_size)


            # self._window = self._window[1:]
            # self._window = self._window.append(self.frame_generator.__next__())





        

                

    def _finalize_run(self) -> Tuple[
        Np2DNumberArray,
    ]:
        """
        Finalizes the MSD analysis after running.
        """

        return (
            self.msd / ( self._number_windows * self.gap )
        )
