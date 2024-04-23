"""
A module containing the diffcalc class. The Diffcalc class is used to calculate the self-diffusion coefficient of given target molecules in given trajectory file(s) with a specific window size and gap. The self-diffusion coefficient is calculated by the mean square displacement (MSD) of the target species and is averaged over all target species. Self-diffusion coefficient is a measure of how fast a particle diffuses in a medium. It is calculated by the following formula:

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
from .exceptions import DiffCalcError, DiffCalcWarning
from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal
from PQAnalysis.core import distance, Cells
from PQAnalysis.traj import Trajectory, check_trajectory_PBC, check_trajectory_vacuum
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.io import TrajectoryReader

class DiffCalc:
    """
    A class for calculating the self-diffusion coefficient as running average of given target molecules in given trajectory file(s) with a specific window size and gap.

    The Diffcalc class is initialized with the provided parameters. The diffcalc analysis can be run by calling the run method. The run method returns the mean-square distances MSD and the component x,y,z of it, the root-mean-square distance RMSD and its components and the self-diffusion coefficient D and the x,y, z components of D.

    The Diffcalc class can be initialized with either a trajectory object or via a TrajectoryReader object. If a trajectory object is given, it is assumed to have a constant topology over all frames! The main difference between the two is that the TrajectoryReader object allows for lazy loading of the trajectory, meaning that the trajectory is only loaded frame by frame when needed. This can be useful for large trajectories that do not fit into memory.
    """

    _use_full_atom_default = False
    _window_size_default = 1000
    _gap_default = 10
    _n_start_default = 0
    _n_stop_default = 10000000000
    _n_frames_default = 0
   

    def __init__(self,
                 traj: Trajectory | TrajectoryReader,
                 target_species: SelectionCompatible,
                 window_size: PositiveInt | None = 500,
                 gap: PositiveInt | None = 1,
                 n_start: PositiveInt | None = None,
                 n_stop: PositiveInt | None = None,
                 n_frames: PositiveInt | None = None,
                 use_full_atom_info: bool = False,
                 ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The trajectory to analyze. If a TrajectoryReader is provided, the trajectory frame by frame via a frame_generator
        target_species : SelectionCompatible
            The target species of the self-diffusion coefficient calculation.
        window_size : PositiveInt, optional
            The window size of the running average, by default 500.
        gap : PositiveInt, optional
            The gap of the running average, by default 1.
        n_start : PositiveInt, optional
            The starting frame of the trajectory, by default None.
        n_stop : PositiveInt, optional
            The stopping frame of the trajectory, by default None.
        n_frames : PositiveInt, optional
            The number of frames to use, by default None.
        use_full_atom_info : bool, optional
            Whether to use the full atom information of the trajectory or not, by default None (False).
    

      

        Raises
        ------
        DiffCalcError
            If the trajectory is not fully periodic or fully in vacuum. Meaning that some frames are in vacuum and others are periodic.
        DiffCalcError
            If the trajectory is empty.
        DiffCalcError
            If target_species is empty.

        Notes
        -----        
        Furthermore, to initialize the diffcalc object target_species has to be defined. The target_species is a selection string that selects the target atoms of the self-diffusion coefficient calculation. The selection string can be any valid selection string for the topology selection. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.

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

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        if window_size is None:
            self.window_size = self._window_size_default
        else:
            self.window_size = window_size
        
        if gap is None:
            self.gap = self._gap_default
        else:
            self.gap = gap
        
        if n_start is None:
            self.n_start = self._n_start_default
        else:
            self.n_start = n_start
        
        if n_stop is None:
            self.n_stop = self._n_stop_default
        else:
            self.n_stop = n_stop
                
        if n_frames is None:
            self.n_frames = self._n_frames_default
        else:
            self.n_frames = n_frames



        ################################
        # Initialize Selection objects #
        ################################

     
        self.target_species = target_species

        if selection =="resid~1":
            residue_names_1 = np.argwhere(np.isin(topology._residue_ids, 1)).flatten()
            residue_numbers_1  = topology._residue_numbers[residue_names_1]
            print(residue_names_1)
            print(residue_numbers_1)
            print(topology._residue_numbers)
            print(topology._residue_ids)

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
            raise DiffCalcError("Trajectory cannot be of length 0.")

        self.first_frame = next(self.frame_generator)
        self.topology = traj.topology

    

        self.target_indices = self.target_selection.select(
            self.topology, self.use_full_atom_info)
        
        self.residue_length = len(self.target_indices)

    def __str__(self) -> str:
        """
        Returns a string representation of the diffcalc object.

        Returns
        -------
        str
            A string representation of the diffcalc object.
        """
        return f"DiffCalc object with window size {self.window_size}, gap {self.gap}, start frame {self.n_start}, stop frame {self.n_stop}, number of configurations {self.n_frames}, target species {self.target_species}, using full atom information {self.use_full_atom_info}"
    
    def __repr__(self) -> str:
        """
        Returns a string representation of the diffcalc object.

        Returns
        -------
        str
            A string representation of the diffcalc object.
        """
        return f"DiffCalc(window_size={self.window_size}, gap={self.gap}, n_start={self.n_start}, n_stop={self.n_stop}, n_frames={self.n_frames}, target_species={self.target_species}, use_full_atom_info={self.use_full_atom_info}"
    def _center_of_mass(self) -> None:
        """
        Calculates the center of mass of the target species.

        This method calculates the center of mass of the target species and writes the center of mass in a new trajectory file.

        """


    def _add_origin(self, active_origins: int):
        """
        Adds an origin to the diffcalc .

        This method adds an origin to the diffcalc and calculates the mean square displacement (MSD) of the target species.

        Parameters
        ----------
        active_origins : int
            The number of active origins.

        """
        self.origins[active_origins] = self.target_pos

    def _remove_origin(self):
        """
        Removes an origin from the diffcalc .

        This method removes an origin from the diffcalc and calculates the mean square displacement (MSD) of the target species.

        """
        self.origins = np.delete(self.origins, 0, axis=0)
        self.image_distances = np.delete(self.image_distances, 0, axis=0)
    def _shift_origin(self):
        """
        Shifts the origin of the diffcalc .

        This method shifts the origin of the diffcalc and calculates the mean square displacement (MSD) of the target species.

        """
        self.origins = np.roll(self.origins, -1, axis=0)
        self.image_distances = np.roll(self.image_distances, -1, axis=0)

    def _check_trajectory_conditions(self):
        """
        Checks if the trajectory is fully periodic or fully in vacuum.

        Raises
        ------
        DiffCalcError
            If the trajectory is not fully periodic or fully in vacuum. Meaning that some frames are in vacuum and others are periodic.
        """
        if not check_trajectory_PBC(self.cells) and not check_trajectory_vacuum(self.cells):
            raise DiffCalcError(
                "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the diffcalc .")
    
    # @timeit_in_class
    def run(self) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
        """
        Runs the diffcalc .

        This method runs the calculation of the diffusion coefficient of the target species and returns the MSD(x,y,z),RMSD(x,y,z),D(x,y,z) of the target species.

        This method will display a progress bar by default. This can be disabled by setting with_progress_bar to False.

        Returns
        -------
   
         """
        self._initialize_run()
        self._calculate_msd()
        return self._finalize_run()
    
    def _initialize_run(self):
        """
        Initializes the diffcalc for running.

        This method is called by the run method of the DiffCalc class. It initializes the number of origins which indicates the number of completing one window and the total number of origins which indicates the total number of completing the used trajectory from starting and end point. The used trajectory is defined by the start_frame and the stop_frame. 
        The method also initializes the targets position, the targets previous position, the image distances and the origins of the diffcalc. Furthermore, the method prints the initialized run parameters.

        """
        if self.n_frames == 0:
            self.n_frames = self.n_stop - self.n_start
      
        
        self.n_origin  = np.fix(self.window_size / self.gap) # number of orgin per window
        self.total_origin = np.fix((self.n_frames - self.window_size - self.n_start) / self.gap) # number of total origin for the used trajectory
        if self.n_start > (self.n_frames - self.window_size):
            raise DiffCalcError("Starting configuration {self.n_start} is larger than the number of configurations {self.n_frames} minus the window size {self.window_size}. \
                                Establish a window size {self.window_size} within the {self.n_frames} total configurations minus the starting configuration {self.n_start}.")
        
        if self.n_start > 0:
            print(f"Starting from configuration {self.n_start} of {self.n_frames} configurations.")




        self.stop_frame = np.floor((self.n_frames - self.window_size) / self.gap) * self.gap # the last frame of the used trajectory


        if self.stop_frame == 0:
            self.stop_frame = 1
        
        if self.topology is None:
            for i,target_index in enumerate(self.target_indices):
                self.target_pos = self.first_frame.pos[target_index]
                self.target_prev_pos = np.zeros((len(self.target_indices), 3))
                self.image_distances = np.zeros((len(self.target_indices), 3))
                self.origins = np.zeros((self.n_origin, len(self.target_indices), 3))
          
        else:
            self._center_of_mass()
            for i,target_index in enumerate(self.target_indices):
                self.target_prev_pos = np.zeros(len(self.target_pos), 3)
                self.image_distances = np.zeros(len(self.target_pos), 3)
                self.origins = np.zeros((self.n_origin, len(self.target_indices), 3))

        print("###### Initialized run ######")
        print(f"Window size: {self.window_size}")
        print(f"Gap: {self.gap}")
        print(f"Start frame: {self.n_start}")
        print(f"Stop frame: {self.n_stop}")
        print(f"Number of configurations: {self.n_frames}")
        print(f"Number of target species: {len(self.target_indices)}")
        print("Cells: ", self.cells)
        print("Target species: ", self.target_species)
        print("Target indices: ", self.target_indices)
     
        print("Using full atom information" if self.use_full_atom_info else "Not using full atom information")
        print("Using topology and center of mass as reference for the target species" if self.topology is not None else "Not using topology as reference the target species are used")




    

    def _calculate_msd(self):
        """
    #     Calculates the mds of the trajectory.
    #     """
        total_frames = 0
        stepping = 0
        active_origins = 0
        last = 0
        origin_counter = 0
        entry = 0
        for frame in tqdm(self.frame_generator, total=self.n_frames, desc="Calculating MSD"):
            if total_frames >= self.n_start:
                stepping = total_frames % self.gap
                if stepping == 0:
                    if (active_origins != self.n_origin) and (total_frames <= self.stop_frame):
                        if active_origins == 0:
                            last = total_frames
                        
                        self._add_origin(active_origins)
                        active_origins += 1
                        origin_counter += 1

                    if (last + self.window_size) == total_frames:
                        # loop over all target species used the distace function to calculate the mean square displacement from this position and last frame position
                        for frame in tqdm(self.frame_generator, total=self.n_frames, disable=not config.with_progress_bar):
                            distances = distance(self.target_pos,self.target_prev_pos,cells=self.cells)
                            distances 
                            self.image_distances = distance(self.image_distances, distances, cells=self.cells)
                            distances = self.target_pos - self.origins + self.image_distances 
                            self.msd = np.sum(self.msd, axis=1) + np.sum(distances**2, axis=1)
                        if total_frames > self.n_stop:
                            self._remove_origin()
                            active_origins -= 1
                            last = last + self.gap
                        else:
                            self._shift_origin()
                            last = last + self.gap
                            origin_counter += 1

                    entry = total_frames - last
                  
                    for i in range(0,active_origins-1,-1):
                        for j in range(0,len(self.target_indices)-1,1):
                            if entry != 0:
                                distances = distance(self.target_pos,self.target_prev_pos,cells=self.cells)

                                self.image_distances = distance(self.image_distances, distances, cells=self.cells)
                            distances = self.target_pos - self.origins + self.image_distances
                            self.msd = np.sum(self.msd, axis=1) + np.sum(distances**2, axis=1)
                        entry = entry - self.gap

        
    def _finalize_run(self):
        """
        Finalizes the diffcalc after running.

        """
        total_origins = self.stop_frame / self.gap
        self.diffusion_coeff = self.msd / len(self.target_indices) / total_origins 
        self.diffusion_coeff = self.diffusion_coeff / 6 / self.gap

        self.diffusion_coeff_x = self.msd_x / len(self.target_indices) / total_origins
        self.diffusion_coeff_x = self.diffusion_coeff_x / 6 / self.gap

        self.diffusion_coeff_y = self.msd_y / len(self.target_indices) / total_origins
        self.diffusion_coeff_y = self.diffusion_coeff_y / 6 / self.gap


    