"""
A module containing the NEPWriter class to 
write NEP training and testing files.
"""

import logging
import _io
import numpy as np

from beartype.typing import List, Dict

from PQAnalysis.io import (
    BaseWriter,
    FileWritingMode,
    OutputFileFormat,
    read_trajectory_generator,
    calculate_frames_of_trajectory_file,
    EnergyFileReader
)
from PQAnalysis.io.virial.api import read_stress_file, read_virial_file
from PQAnalysis.utils.units import kcal_per_mol, eV
from PQAnalysis.utils.files import find_files_with_prefix
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.utils.random import get_random_seed
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj import Trajectory, TrajectoryFormat
from PQAnalysis.types import PositiveReal
from PQAnalysis import config
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import NEPError



class NEPWriter(BaseWriter):

    """
    A class to write NEP training and testing files.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None,
        mode: FileWritingMode | str = "w",
    ) -> None:
        """
        Parameters
        ----------
        filename : str
            _description_
        mode : FileWritingMode | str, optional
            _description_, by default "w"
        """
        super().__init__(filename, mode)

        # set original_mode - baseWriter modifies mode after opening the file
        self.original_mode = FileWritingMode(mode)

        #########################
        # dummy initializations #
        #########################

        self.use_forces = False
        self.use_stress = False
        self.use_virial = False

        self.xyz_file_extension = None
        self.energy_file_extension = None
        self.info_file_extension = None
        self.force_file_extension = None
        self.stress_file_extension = None
        self.virial_file_extension = None

        self.train_file = None
        self.train_writer = None
        self.test_file = None
        self.test_writer = None
        self.validation_file = None
        self.validation_writer = None
        self.validation_ref_file = None
        self.validation_ref_writer = None

        self.n_train_frames = 0
        self.n_test_frames = 0
        self.n_validation_frames = 0

        self.test_ratio = 0.0
        self.validation_ratio = 0.0
        self.is_validation = False

    @runtime_type_checking
    def write_from_files(
        self,
        file_prefixes: List[str] | str,
        use_forces: bool = False,
        use_stress: bool = False,
        use_virial: bool = False,
        xyz_file_extension: str | None = None,
        energy_file_extension: str | None = None,
        info_file_extension: str | None = None,
        force_file_extension: str | None = None,
        stress_file_extension: str | None = None,
        virial_file_extension: str | None = None,
        test_ratio: PositiveReal = 0.0,
        total_ratios: str | None = None,
    ) -> None:
        """
        Writes the NEP trajectory file from the given files.

        Parameters
        ----------
        file_prefixes : List[str] | str
            The prefixes of the files to find. Here with prefix we mean
            the part of the filename not only the name before the
            extension, but every matching file that starts with the
            given prefix.
        use_forces : bool, optional
            Whether to include forces in the output file, by default False
        use_stress : bool, optional
            Whether to include the stress tensor in the output file,
            by default False
        use_virial : bool, optional
            Whether to include the virial in the output file, by default False
        xyz_file_extension : str, optional
            The extension of the xyz files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with the
            given file prefixes.
        energy_file_extension : str, optional
            The extension of the energy files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with the
            given file prefixes.
        info_file_extension : str, optional
            The extension of the info files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with
            the given file prefixes.
        force_file_extension : str, optional
            The extension of the force files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with the
            given file prefixes.
        stress_file_extension : str, optional
            The extension of the stress files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with the
            given file prefixes.
        virial_file_extension : str, optional
            The extension of the virial files, by default None.
            This means that the respective file extension will
            be automatically determined from all files with the
            given file prefixes.
        test_ratio : PositiveReal, optional
            The ratio of testing frames to the total number of
            frames, by default 0.0. If the test_ratio is 0.0 no
            train and test files are created. If the test_ratio
            is larger not equal to 0.0, the test_ratio is used
            to determine the number of training and testing
            frames. The final ratio will be as close to the
            test_ratio as possible, but if it is not possible
            to have the exact ratio, always the higher next
            higher ratio is chosen. As output filenames the
            original filename is used with the suffix _train
            or _test appended and the same FileWritingMode as
            the original file is used.
        total_ratios: str, optional
            The total_ratios keyword argument is used to describe
            frame ratios including validation frames in the format
            train_ratio:test_ratio:validation_ratio. The
            validation_ratio is optional and if not given, no
            validation frames are written. The total sum of the
            integer values provided do not have to add up to the
            total number of frames in the input trajectory files.
            The ratios are used to determine the ratios of the
            training, testing, and validation frames. The final
            ratio will be as close to the given ratios as possible,
            but if it is not possible to have the exact ratio,
            always the next higher ratio is chosen. As output
            filenames the original filename is used with the suffix
            _train, _test, or _validation appended and the same
            FileWritingMode as the original file is used.
            The validation frames are written to a file with the
            suffix _validation and a file with the suffix
            _validation.ref. The _validation file contains only
            the coordinates and box information to function as
            crude testing input and the _validation.ref file
            contains all information additionally provided in
            the original files.
            Pay Attention: This keyword argument is mutually
            exclusive with the test_ratio keyword argument.
            If both are given, a ValueError is raised.

        Raises
        ------
        ValueError
            If the test_ratio is larger than 1.0.
        ValueError
            If test_ratio and total_ratios are given at the same time.
        """

        self.use_forces = use_forces
        self.use_stress = use_stress
        self.use_virial = use_virial

        self.xyz_file_extension = xyz_file_extension
        self.energy_file_extension = energy_file_extension
        self.info_file_extension = info_file_extension
        self.force_file_extension = force_file_extension
        self.stress_file_extension = stress_file_extension
        self.virial_file_extension = virial_file_extension

        files = find_files_with_prefix(file_prefixes)

        file_dict = self._determine_files(
            files,
            file_prefixes,
        )

        xyz_files = file_dict[OutputFileFormat.XYZ.value]
        en_files = file_dict[OutputFileFormat.ENERGY.value]
        info_files = file_dict[OutputFileFormat.INFO.value]
        force_files = file_dict[OutputFileFormat.FORCE.value
                                ] if use_forces else []
        stress_files = file_dict[OutputFileFormat.STRESS.value
                                 ] if use_stress else []
        virial_files = file_dict[OutputFileFormat.VIRIAL.value
                                 ] if use_virial else []

        if not np.isclose(test_ratio, 0.0) or total_ratios is not None:
            if self.filename is None:
                self.logger.error(
                    (
                        "No output filename was specified. In order to "
                        "use the test_ratio or total_ratios splitting feature "
                        "a filename has to be given in order to write all "
                        "generated files."
                    ),
                    exception=NEPError
                )

            self._setup_frame_splitting_for_training(test_ratio, total_ratios)

        self.open()

        self.n_train_frames = 0
        self.n_test_frames = 0
        self.n_validation_frames = 0

        for i, xyz_file in enumerate(xyz_files):

            stress = read_stress_file(stress_files[i]) if use_stress else []
            virial = read_virial_file(virial_files[i]) if use_virial else []
            energy = EnergyFileReader(en_files[i], info_files[i]).read()

            xyz_generator = read_trajectory_generator(
                xyz_file, traj_format=TrajectoryFormat.XYZ
            )
            force_generator = read_trajectory_generator(
                force_files[i], traj_format=TrajectoryFormat.FORCE
            ) if use_forces else None

            n_frames = calculate_frames_of_trajectory_file(xyz_file)

            n_train_max, n_test_max, _ = self._get_effective_training_portions(
                n_frames,
                self.test_ratio,
                self.validation_ratio,
            )

            seed = get_random_seed()
            rng = np.random.default_rng(seed=seed)

            indices = np.arange(n_frames)
            rng.shuffle(indices)

            _train = n_train_max - self.n_train_frames
            _test = n_test_max - self.n_test_frames

            train_indices = indices[:_train]
            test_indices = indices[_train:_train + _test]
            validation_indices = indices[_train + _test:]

            for j, system in enumerate(xyz_generator):
                if hasattr(energy, "qm_energy"):
                    system.energy = energy.qm_energy[j]
                else:
                    self.logger.error(
                        (
                            "No QM energy found in the energy file. "
                            "The NEP builder is implemented only for QM energies. "
                            "If there is a need for a different energy type, "
                            "please file an issue at "
                            f"{config.code_base_url}issues."
                        ),
                        exception=NEPError
                    )

                if use_forces:
                    force_system = next(force_generator)
                    system.forces = force_system.forces

                if use_virial and virial:
                    system.virial = virial[j]
                if use_stress and stress:
                    system.stress = stress[j]

                self.write_from_atomic_system(
                    system, self.file, use_forces, use_stress, use_virial
                )

                self.is_validation = False

                file_to_write = None  # to avoid linting warning

                if j in train_indices:
                    self.n_train_frames += 1
                    file_to_write = self.train_file
                elif j in test_indices:
                    self.n_test_frames += 1
                    file_to_write = self.test_file
                elif j in validation_indices:
                    self.n_validation_frames += 1
                    self.is_validation = True
                    file_to_write = self.validation_ref_file
                else:
                    self.logger.error(
                        (
                            "An error occurred during the frame splitting. "
                            "The frame index is not in any of the training, "
                            "testing, or validation indices. This should not "
                            "happen. Please file an issue at "
                            f"{config.code_base_url}issues."
                        ),
                        exception=NEPError
                    )

                self.write_from_atomic_system(
                    system,
                    file_to_write,
                    use_forces,
                    use_stress,
                    use_virial,
                )

                if self.is_validation:
                    self.write_from_atomic_system(
                        system,
                        self.validation_file,
                        use_forces=False,
                        use_stress=False,
                        use_virial=False,
                    )

            self.logger.info(
                (
                    f"Processed {n_frames} frames from files:\n"
                    f"{xyz_file}, "
                    f"{en_files[i]}, "
                    f"{info_files[i]}, "
                    f"{force_files[i] if use_forces else None}, "
                    f"{stress_files[i] if use_stress else None}, "
                    f"{virial_files[i] if use_virial else None}"
                )
            )

        self._close_files()

    def _close_files(self):
        """
        Closes the files used for writing the NEP trajectory file.
        """

        if self.train_file is not None:
            self.train_writer.close()
            self.test_writer.close()

        if self.validation_file is not None:
            self.validation_writer.close()
            self.validation_ref_writer.close()

        self.close()

    def _setup_frame_splitting_for_training(
        self,
        test_ratio: PositiveReal = 0.0,
        total_ratios: str | None = None,
    ) -> None:
        """
        Sets up the frame splitting for training.

        It determines the number of training, testing, and validation frames
        based on the given test_ratio or total_ratios keyword arguments.
        It also sets up the files to write the training, testing, and
        validation frames to.

        Parameters
        ----------
        test_ratio : PositiveReal, optional
            The ratio of testing frames to the total number of frames,
            by default 0.0
        total_ratios : str | None, optional
            The total_ratios keyword argument is used to describe
            frame ratios including validation frames in the format
            train_ratio:test_ratio:validation_ratio. The validation_ratio
            is optional and if not given, no validation frames are
            written. The total sum of the integer values provided do
            not have to add up to the total number of frames in the
            input trajectory files. The ratios are used to determine
            the ratios of the training, testing, and validation
            frames. The final ratio will be as close to the given
            ratios as possible, but if it is not possible to have
            the exact ratio, always the next higher ratio is chosen.
            As output filenames the original filename is used with
            the suffix _train, _test, or _validation appended and the
            same FileWritingMode as the original file is used.


        Raises
        ------
        ValueError
            If the test_ratio and total_ratios keyword arguments
            are mutually exclusive.
        ValueError
            If the total_ratios keyword argument is not in the correct format.
        ValueError
            If the test_ratio is between 0.0 and 1.0.
        ValueError
            If validation frames are given without test frames.
        """

        if not np.isclose(test_ratio, 0.0) and total_ratios is not None:
            self.logger.error(
                (
                    "The test_ratio and total_ratios keyword "
                    "arguments are mutually exclusive."
                ),
                exception=NEPError
            )

        n_train = 0.0
        n_validation = 0.0
        n_test = 0.0

        if total_ratios is not None:
            ratios = total_ratios.split(":")
            if len(ratios) == 2:
                n_train = float(ratios[0])
                n_test = float(ratios[1])
                n_validation = 0.0
            elif len(ratios) == 3:
                n_train = float(ratios[0])
                n_test = float(ratios[1])
                n_validation = float(ratios[2])
            else:
                self.logger.error(
                    (
                        f"The total_ratios keyword argument {total_ratios} "
                        "is not in the correct format. The correct format "
                        "is train_ratio:test_ratio:validation_ratio."
                    ),
                    exception=NEPError
                )

        sum_frames = n_train + n_test + n_validation

        self.test_ratio = n_test / sum_frames
        self.validation_ratio = n_validation / sum_frames

        if self.test_ratio > 1.0:
            self.logger.error(
                (
                    "The test_ratio must be between 0.0 and 1.0. "
                    "The given test_ratio "
                    f"is {self.test_ratio}."
                ),
                exception=NEPError
            )

        if np.isclose(self.test_ratio,
                      0.0) and not np.isclose(self.validation_ratio, 0.0):
            self.logger.error(
                (
                    "It has no sense to have validation frames without test "
                    "frames. This error results from the given total_ratios "
                    f"keyword argument {total_ratios}."
                ),
                exception=NEPError
            )

        self.train_file = None
        self.train_writer = None
        self.train_file = None
        self.train_writer = None
        self.validation_file = None
        self.validation_writer = None
        self.validation_ref_file = None
        self.validation_ref_writer = None

        if not np.isclose(self.test_ratio, 0.0):
            self.train_writer = BaseWriter(
                f"{self.filename}_train", mode=self.original_mode
            )
            self.test_writer = BaseWriter(
                f"{self.filename}_test", mode=self.original_mode
            )
            self.train_writer.open()
            self.test_writer.open()
            self.train_file = self.train_writer.file
            self.test_file = self.test_writer.file

        if not np.isclose(self.validation_ratio, 0.0):
            self.validation_writer = BaseWriter(
                f"{self.filename}_validation", mode=self.original_mode
            )
            self.validation_writer.open()
            self.validation_ref_writer = BaseWriter(
                f"{self.filename}_validation.ref", mode=self.original_mode
            )
            self.validation_ref_writer.open()
            self.validation_file = self.validation_writer.file
            self.validation_ref_file = self.validation_ref_writer.file

    def _determine_files(
        self,
        files: List[str],
        file_prefixes: List[str],
    ) -> Dict[str, List[str]]:
        """
        Determines the files to be used for writing the NEP trajectory file.

        Parameters
        ----------
        files : List[str]
            The files to be used for writing the NEP trajectory file.
        file_prefixes : List[str]
            The prefixes of the files to find.

        Returns
        -------
        Dict[str, List[str]]
            The files to be used for writing the NEP trajectory file.
            The keys are the file types and the values 
            are the respective files.

        Raises
        ------
        ValueError
            If the number of files does not match the 
            number of coordinate files.
        """

        file_dict = {}

        file_dict[OutputFileFormat.XYZ.value] = self._get_files(
            files,
            OutputFileFormat.XYZ,
            self.xyz_file_extension,
            file_prefixes
        )

        file_dict[OutputFileFormat.ENERGY.value] = self._get_files(
            files,
            OutputFileFormat.INSTANTANEOUS_ENERGY,
            self.energy_file_extension,
            file_prefixes
        )

        file_dict[OutputFileFormat.INFO.value] = self._get_files(
            files,
            OutputFileFormat.INFO,
            self.info_file_extension,
            file_prefixes
        )

        if self.use_forces:
            file_dict[OutputFileFormat.FORCE.value] = self._get_files(
                files,
                OutputFileFormat.FORCE,
                self.force_file_extension,
                file_prefixes
            )

        if self.use_stress:
            file_dict[OutputFileFormat.STRESS.value] = self._get_files(
                files,
                OutputFileFormat.STRESS,
                self.stress_file_extension,
                file_prefixes
            )

        if self.use_virial:
            file_dict[OutputFileFormat.VIRIAL.value] = self._get_files(
                files,
                OutputFileFormat.VIRIAL,
                self.virial_file_extension,
                file_prefixes
            )

        def raise_number_of_files_error(
            file_type: str,
            files: List[str],
            xyz_files: List[str],
        ):
            """
            Raises an error if the number of files does not 
            match the number of coordinate files.

            Parameters
            ----------
            file_type : str
                the type of the files
            files : List[str]
                the files to be used for writing the NEP trajectory file
            xyz_files : List[str]
                the coordinate files
            """
            self.logger.error(
                (
                    f"The number of {file_type} files does not match "
                    "the number of coordinate files. The found "
                    f"{file_type} files are: {files} "
                    "and the found coordinate files are: "
                    f"{xyz_files}"
                ),
                exception=NEPError
            )

        en_files = file_dict[OutputFileFormat.ENERGY.value]
        xyz_files = file_dict[OutputFileFormat.XYZ.value]
        info_files = file_dict[OutputFileFormat.INFO.value]
        force_files = file_dict[OutputFileFormat.FORCE.value
                                ] if self.use_forces else []
        stress_files = file_dict[OutputFileFormat.STRESS.value
                                 ] if self.use_stress else []
        virial_files = file_dict[OutputFileFormat.VIRIAL.value
                                 ] if self.use_virial else []

        if len(en_files) != len(xyz_files):
            raise_number_of_files_error("energy", en_files, xyz_files)

        if len(info_files) != len(xyz_files):
            raise_number_of_files_error("info", info_files, xyz_files)

        if self.use_forces and len(force_files) != len(xyz_files):
            raise_number_of_files_error("force", force_files, xyz_files)

        if self.use_stress and len(stress_files) != len(xyz_files):
            raise_number_of_files_error("stress", stress_files, xyz_files)

        if self.use_virial and len(virial_files) != len(xyz_files):
            raise_number_of_files_error("virial", virial_files, xyz_files)

        self.logger.info(
            f"""
Reading files to write NEP trajectory file:
    - xyz_files:    {xyz_files}
    - en_files:     {en_files}
    - info_files:   {info_files}
    - force_files:  {force_files if self.use_forces else ""}
    - stress_files: {stress_files if self.use_stress else ""}
    - virial_files: {virial_files if self.use_virial else ""}
"""
        )

        return file_dict

    def _get_files(
        self,
        files: List[str],
        output_file_format: OutputFileFormat,
        file_extension: str | None,
        file_prefixes: List[str],
    ) -> List[str]:
        """
        Get the files that match the given file prefixes and file extension.

        Parameters
        ----------
        files : List[str]
            The files to be used for writing the NEP trajectory file.
        outputFileFormat : OutputFileFormat
            The file format of the files to be used for 
            writing the NEP trajectory file.
        file_extension : str | None
            The file extension of the files to be used for
            writing the NEP trajectory file.
        file_prefixes : List[str]
            The prefixes of the files to find.

        Returns
        -------
        List[str]
            The files that match the given file prefixes and file extension.

        Raises
        ------
        ValueError
            If no files with the given file prefixes were found.
        ValueError
            If no files with the given file prefixes 
            and file extension were found.
        """

        filtered_files = OutputFileFormat.find_matching_files(
            files, output_file_format, file_extension
        )

        if len(filtered_files) == 0:
            if file_extension is not None:
                self.logger.error(
                    (
                        "You did specify a file extension for the "
                        f"{output_file_format} files, but no files with "
                        f"the extension {file_extension} were found, "
                        f"that match the given file prefixes {file_prefixes}."
                    ),
                    exception=NEPError
                )
            else:
                output_file_formats = OutputFileFormat.get_file_extensions(
                    output_file_format
                )
                self.logger.error(
                    (
                        f"No {output_file_format} files were found in "
                        f"{files} that match the given file prefixes "
                        f"{file_prefixes}. All possible file extensions are "
                        f"{output_file_formats}. "
                        "If the specific file extension you are looking for "
                        "is not in the list, please specify it using the "
                        "corresponding file_extension argument. If the "
                        "files should be found, please check the file paths "
                        "and the file prefixes. Additionally, if you think "
                        "that the file extension you chose is of general "
                        "interest and should be added to the list of "
                        "possible file extensions, please file an issue at "
                        f"{config.code_base_url}issues."
                    ),
                    exception=NEPError
                )

        return sorted(filtered_files)

    def _get_effective_training_portions(
        self,
        n_frames: int,
        test_ratio: PositiveReal,
        validation_ratio: PositiveReal,
    ) -> tuple[int, int, int]:
        """
        Calculates the maximum number of training, testing,
        and validation frames.

        By calculating the maximum number of training, testing,
        and validation frames, the number of frames to be added 
        to the training and testing files can be determined.


        Parameters
        ----------
        n_frames : int
            the number of frames to add to the training and testing files
        test_ratio : PositiveReal
            the ratio of testing frames to the total number of frames
        validation_ratio : PositiveReal
            the ratio of validation frames to the total number of frames

        Returns
        -------
        tuple[int, int, int]
            the maximum number of training and testing frames
        """

        detected_frames = self.n_train_frames + \
            self.n_test_frames + self.n_validation_frames
        total_frames = n_frames + detected_frames

        max_test_frames = int(np.ceil(total_frames * test_ratio))
        max_validation_frames = int(np.ceil(total_frames * validation_ratio))
        max_train_frames = total_frames - max_test_frames - max_validation_frames

        return max_train_frames, max_test_frames, max_validation_frames

    @runtime_type_checking
    def write_from_trajectory(
        self,
        trajectory: Trajectory,
        use_forces: bool = False,
        use_stress: bool = False,
        use_virial: bool = False,
    ) -> None:
        """
        Writes the NEP trajectory file from the given trajectory.

        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to be written to the NEP trajectory file.
        use_forces : bool, optional
            Whether to write the forces to the NEP trajectory file, 
            by default False
        use_stress : bool, optional
            Whether to write the stress tensor to the NEP trajectory file, 
            by default False
        use_virial : bool, optional
            Whether to write the virial tensor to the NEP trajectory file, 
            by default False
        """

        self.open()
        for frame in trajectory:
            self.write_from_atomic_system(
                frame, use_forces, use_stress, use_virial
            )

        self.close()

    @runtime_type_checking
    def write_from_atomic_system(
        self,
        system: AtomicSystem,
        file: _io.TextIOWrapper | None,
        use_forces: bool = False,
        use_stress: bool = False,
        use_virial: bool = False,
    ) -> None:
        """
        Writes the NEP trajectory file from the given atomic system.

        Parameters
        ----------
        system : AtomicSystem
            The system to be written to the NEP trajectory file.
        file : _io.TextIOWrapper
            The file to write the NEP trajectory file to.
            If None, nothing is written
        use_forces : bool, optional
            Whether to write the forces to the NEP trajectory file, 
            by default False
        use_stress : bool, optional
            Whether to write the stress tensor to the NEP trajectory file,
            by default False
        use_virial : bool, optional
            Whether to write the virial tensor to the NEP trajectory file,
            by default False

        Raises
        ------
        ValueError
            If the system does not have coordinates.
        ValueError
            If the system does not have an energy.
        ValueError
            If the system does not have forces and they
            were specified to be written to the NEP trajectory file.
        ValueError
            If both the stress and the virial tensor were
            specified to be written to the NEP trajectory file.
            Only one of them can be written at a time.
        ValueError
            If the system does not have a stress tensor and it
            was specified to be written to the NEP trajectory file.
        ValueError
            If the system does not have a virial tensor and it
            was specified to be written to the NEP trajectory file.
        """

        if file is None:
            return

        if not system.has_pos:
            self.logger.error(
                (
                    "The system does not have coordinates, "
                    "which are required for NEP trajectory files."
                ),
                exception=NEPError
            )

        if not system.has_energy:
            self.logger.error(
                "The system does not have an energy, "
                "which is required for NEP trajectory files.",
                exception=NEPError
            )

        if use_forces and not system.has_forces:
            self.logger.error(
                (
                    "The system does not have forces, and they were "
                    "specified to be written to the NEP trajectory file."
                ),
                exception=NEPError
            )

        if use_stress and use_virial:
            self.logger.error(
                (
                    "Both the stress and the virial tensor were "
                    "specified to be written to the NEP trajectory file. "
                    "Only one of them can be written at a time."
                ),
                exception=NEPError
            )

        if use_stress and not system.has_stress:
            self.logger.error(
                (
                    "The system does not have a stress tensor, "
                    "and it was specified to be written to the NEP trajectory file."
                ),
                exception=NEPError
            )

        if use_virial and not system.has_virial:
            self.logger.error(
                (
                    "The system does not have a virial tensor, and "
                    "it was specified to be written to the NEP trajectory file."
                ),
                exception=NEPError
            )

        self._write_header(system, file, use_forces, use_stress, use_virial)
        self._write_body(system, file, use_forces)

    def _write_header(
        self,
        system: AtomicSystem,
        file: _io.TextIOWrapper,
        use_forces: bool = False,
        use_stress: bool = False,
        use_virial: bool = False,
    ) -> None:
        """
        Writes the header of the NEP trajectory file.

        Parameters
        ----------
        system : AtomicSystem
            The system to be written to the NEP trajectory file.
        file : _io.TextIOWrapper
            The file to write the NEP trajectory file to.
        use_forces : bool, optional
            Whether to write the forces to the NEP trajectory file,
            by default False
        use_stress : bool, optional
            Whether to write the stress tensor to the NEP trajectory file,
            by default False
        use_virial : bool, optional
            Whether to write the virial tensor to the NEP trajectory file, 
            by default False
        """

        box_matrix = np.transpose(system.cell.box_matrix)

        print(system.n_atoms, file=file)

        energy = system.energy * eV / kcal_per_mol

        file.write(f"energy={energy} ")

        file.write("config_type=nep2xyz ")

        file.write("lattice=\"")
        for i in range(3):
            for j in range(3):
                file.write(f"{box_matrix[i][j]} ")
        file.write("\" ")

        if use_virial:
            virial = system.virial * eV / kcal_per_mol

            file.write("virial=\"")
            for i in range(3):
                for j in range(3):
                    file.write(f"{virial[i][j]} ")
            file.write("\" ")

        if use_stress:
            stress = system.stress * eV / kcal_per_mol

            file.write("stress=\"")
            for i in range(3):
                for j in range(3):
                    file.write(f"{stress[i][j]} ")
            file.write("\" ")

        file.write("properties=species:S:1:pos:R:3")
        if use_forces:
            file.write(":forces:R:3")
        file.write("\n")

    def _write_body(
        self,
        system: AtomicSystem,
        file: _io.TextIOWrapper,
        use_forces: bool = False,
    ) -> None:
        """
        Writes the body of the NEP trajectory file.

        Parameters
        ----------
        system: AtomicSystem
            The system to be written to the NEP trajectory file.
        file : _io.TextIOWrapper
            The file to write the NEP trajectory file to.
        use_forces: bool, optional
            Whether to write the forces to the NEP trajectory file, by default False
        """

        for i in range(system.n_atoms):
            symbol = system.atoms[i].symbol
            symbol = symbol[0].upper() + symbol[1:]

            print(
                (
                    f"{symbol:<4} "
                    f"{system.pos[i][0]:12.8f} "
                    f"{system.pos[i][1]:12.8f} "
                    f"{system.pos[i][2]:12.8f}"
                ),
                file=file,
                end=" "
            )

            forces = system.forces * eV / kcal_per_mol

            if use_forces:
                print(
                    (
                        f"{forces[i][0]:15.8e} "
                        f"{forces[i][1]:15.8e} "
                        f"{forces[i][2]:15.8e}"
                    ),
                    file=file,
                    end=" "
                )

            print(file=file)
