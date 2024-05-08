"""
A module containing the MSD class. The MSD class is used to calculate the self-diffusion coefficient of given target molecules in given trajectory file(s) with a specific window size and gap. The self-diffusion coefficient is calculated by the mean square displacement (MSD) of the target species and is averaged over all target species. Self-diffusion coefficient is a measure of how fast a particle diffuses in a medium. It is calculated by the following formula:

.. math::
    D = \\frac{1}{6t} \\langle \\sum_{i=1}^{N} (r_i(t) - r_i(0))^2 \\rangle

where :math:`D` is the self-diffusion coefficient, :math:`t` is the time, :math:`N` is the number of target species, :math:`r_i(t)` is the position of the target species at time :math:`t` and :math:`r_i(0)` is the position of the target species at time :math:`0`.



"""

from __future__ import annotations

# 3rd party imports
import numpy as np
import warnings

# 3rd party imports
from beartype.typing import Tuple, List
from tqdm.auto import tqdm

# local imports
import PQAnalysis.config as config

# local imports
from .exceptions import MSDError, MSDWarning
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cells
from PQAnalysis.traj import Trajectory, check_trajectory_PBC, check_trajectory_vacuum
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.io import TrajectoryReader

class MSD:
    """
    A class for calculating the self-diffusion coefficient as running average of given target molecules in given trajectory file(s) with a specific window size and gap.

    The MSD class is initialized with the provided parameters. The MSD analysis can be run by calling the run method. The run method returns the mean-square distances MSD and the component x,y,z of it, the root-mean-square distance RMSD and its components and the self-diffusion coefficient D and the x,y, z components of D.

    The MSD class can be initialized with either a trajectory object or via a TrajectoryReader object. If a trajectory object is given, it is assumed to have a constant topology over all frames! The main difference between the two is that the TrajectoryReader object allows for lazy loading of the trajectory, meaning that the trajectory is only loaded frame by frame when needed. This can be useful for large trajectories that do not fit into memory.
    """

    _use_full_atom_default = False
    _window_size_default = 1000
    _gap_default = 10
    _frame_start_default = 0
    
   

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
            The trajectory to analyze. If a TrajectoryReader is provided, the trajectory frame by frame via a frame_generator
        target_species : SelectionCompatible
            The target species of the self-diffusion coefficient calculation.
        use_full_atom_info : bool, optional
            Whether to use the full atom information of the trajectory or not, by default None (False).
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
            If the trajectory is not fully periodic or fully in vacuum. Meaning that some frames are in vacuum and others are periodic.
        MSDError
            If the trajectory is empty.
        MSDError
            If target_species is empty.

        Notes
        -----        
        Furthermore, to initialize the MSD object target_species has to be defined. The target_species is a selection string that selects the target atoms of the self-diffusion coefficient calculation. The selection string can be any valid selection string for the topology selection. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.

        See Also
        --------
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        :py:class:`~PQAnalysis.topology.selection.Selection`
        :py:class:`~PQAnalysis.io.trajectoryReader.TrajectoryReader`
        :py:class:`~PQAnalysis.traj.trajectory.Trajectory`
        """
        #####################################################
        # Initialize parameters with default values if None #
        #####################################################        

   
        if window_size is None:
            self.window_size = self._window_size_default
        else:
            self.window_size = window_size
        
        if gap is None:
            self.gap = self._gap_default
        else:
            self.gap = gap

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info
        if no_intra_molecular is None:
            self.no_intra_molecular = self._no_intra_molecular_default
        else:
            self.no_intra_molecular = no_intra_molecular

        if frame_start is None:
            self.frame_start = self._frame_start_default
        else:
            self.frame_start = frame_start

        if n_frames is None and frame_stop is None:
            raise MSDError(
                "Either n_frames or frame_stop has to be defined.")
        elif n_frames is not None and frame_stop is not None:
            raise MSDError(
                "Either n_frames or frame_stop has to be defined, not both.")
        elif n_frames is not None:
            self.n_frames = n_frames
            self.frame_stop = None
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
            # lazy loading of trajectory from file(s)
            self.frame_generator = traj.frame_generator()
        elif len(traj) > 0:
            # use trajectory object as iterator
            self.frame_generator = iter(traj)
        else:
            raise MSDError("Trajectory cannot be of length 0.")

        self.first_frame = next(self.frame_generator)
        self.topology = traj.topology


        self.target_indices = self.target_selection.select(
            self.topology, self.use_full_atom_info)
        
        self.residue_length = len(self.target_indices)

    def __str__(self) -> str:
        """
        Returns a string representation of the MSD object.

        Returns
        -------
        str
            A string representation of the MSD object.
        """
        return f"MSD object with window size {self.window_size}, gap {self.gap}, start frame {self.frame_start}, stop frame {self.frame_stop}, number of configurations {self.n_frames}, target species {self.target_species}, using full atom information {self.use_full_atom_info}"
    
    def __repr__(self) -> str:
        """
        Returns a string representation of the MSD object.

        Returns
        -------
        str
            A string representation of the MSD object.
        """
        return f"MSD(window_size={self.window_size}, gap={self.gap}, frame_start={self.frame_start}, frame_stop={self.frame_stop}, n_frames={self.n_frames}, target_species={self.target_species}, use_full_atom_info={self.use_full_atom_info}"
    



    def _check_trajectory_conditions(self):
        """
        Checks if the trajectory is fully periodic or fully in vacuum.

        Raises
        ------
        MSDError
            If the trajectory is not fully periodic or fully in vacuum. Meaning that some frames are in vacuum and others are periodic.
        """
        if not check_trajectory_PBC(self.cells) and not check_trajectory_vacuum(self.cells):
            raise MSDError(
                "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the MSD .")
    
    # @timeit_in_class
    def run(self) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
        """
        Runs the MSD .

        This method runs the calculation of the diffusion coefficient of the target species and returns the MSD(x,y,z),RMSD(x,y,z),D(x,y,z) of the target species.

        This method will display a progress bar by default. This can be disabled by setting with_progress_bar to False.

        Returns
        -------
   
         """
        self._initialize_run()
        self._calculate_msd()
        return self._finalize_run()
    
    def _initialize_run(self) -> None:
        """
        Initializes the MSD run.

        This method initializes the MSD run by setting the initial values for the MSD calculations.

        """
        self._check_trajectory_conditions()
        self._initialize_arrays()

    def _initialize_arrays(self) -> None:
        """
        Initializes the arrays for the MSD calculations.

        This method initializes the arrays for the MSD calculations.

        """
        
        self.msd_x = np.zeros(self.residue_length)
        self.msd_y = np.zeros(self.residue_length)
        self.msd_z = np.zeros(self.residue_length)
        self.origin_x = np.zeros(self.residue_length)
        self.origin_y = np.zeros(self.residue_length)
        self.origin_z = np.zeros(self.residue_length)
        self.im_x = np.zeros(self.residue_length)
        self.im_y = np.zeros(self.residue_length)
        self.im_z = np.zeros(self.residue_length)

   
    def _calculate_msd(self) -> None:
        """
        Calculates the MSD of the target species by using running average. The window size is defined by the window_size parameter and the gap is defined by the gap parameter. The window is the number of frames to average over and the gap specifies the spacing between the frames in the window.  

        If the the trajectory is fully periodic, the distance between the target species is calculated by the minimum image convention. If the trajectory is fully in vacuum, the distance between the target species is calculated by the Euclidean distance.

        This method calculates the MSD of the target species.

        """
        prev_index = 0
        step_index = 0
        n_steps = np.floor(self.window_size / self.gap)
        n_windows = np.floor(self.n_frames / self.window_size)
        n_total = n_windows * self.window_size
        present_index = 0

        if self.frame_start + self.window_size > self.n_frames:
            raise MSDError("The window size is too large for the trajectory.")
        self
        for frame_index, frame in enumerate(tqdm(self.frame_generator, total=self.n_frames)):
            if frame_index < self.frame_start:
                continue
            if self.frame_stop is not None and frame_index >= self.frame_stop - 1:
                break
            
            if ((frame_index + 1) % self.gap == 0):
                if (step_index != n_steps) and (frame_index != self.frame_stop):
                    if step_index == 0:
                        present_index = frame_index
                    previous_frame = frame
                    step_index += 1
                    n_prev += 1
                if (present_index + self.window_size) == frame_index:
                    present_index = frame_index
                    step_index = 0
                    n_prev = 0
                    continue
                if step_index == n_steps:
                    present_index = frame_index
                    step_index = 0
                    n_prev = 0
                    continue
                if step_index == 0:
                    continue

                if check_trajectory_PBC(self.cells):
                    distance_x, distance_y, distance_z = distance(
                        previous_frame, frame, self.cells)
                elif check_trajectory_vacuum(self.cells):
                    distance_x, distance_y, distance_z = distance(
                        previous_frame, frame, None)
                else:
                    raise MSDError(
                        "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the MSD .")
                self.msd_x += distance_x
                self.msd_y += distance_y
                self.msd_z += distance_z
                self.origin_x += distance_x
                self.origin_y += distance_y
                self.origin_z += distance_z
                self.im_x += distance_x
                self.im_y += distance_y
                self.im_z += distance_z
                prev_index = frame_index
                step_index = 0
                n_prev = 0
                continue
        

        if step_index != 0:
            warnings.warn(
                "The last window is not complete and will be ignored.", MSDWarning)
            
    def _finalize_run(self) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
        """
        Finalizes the MSD run.

        This method finalizes the MSD run by calculating the MSD of the target species.

        Returns
        -------
        Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray,]
            The MSD(x,y,z) of the target species.
        """
        # self.msd_x /= n_total
        # self.msd_y /= n_total
        # self.msd_z /= n_total
        return self.msd_x, self.msd_y, self.msd_z




