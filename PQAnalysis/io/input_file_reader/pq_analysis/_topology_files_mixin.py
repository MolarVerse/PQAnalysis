"""
A mixin class for the input file readers of the analyses that need the
restart file and the moldescriptor file (the topology files) for some keys.
"""

from pathlib import Path

from PQAnalysis.io.input_file_reader.exceptions import InputFileError

from ._parse import _parse_files



class _TopologyFilesMixin:

    """
    A mixin class that infers the topology files (restart file and
    moldescriptor file) from the first trajectory file, if a key requires
    them and they are not given. It needs the attributes and methods of the
    PQAnalysisInputFileReader and a class logger.
    """

    _moldescriptor_file_default = "moldescriptor.dat"

    def _require_topology_files(self, key: str) -> None:
        """
        Requires both topology files for a key, inferring missing ones.

        Parameters
        ----------
        key : str
            The key (or keys) that require the topology files, named in
            the error message.

        Raises
        ------
        InputFileError
            If the files are not given and cannot be inferred.
        """

        if not self._infer_required_topology_files():
            self.logger.error(
                self._missing_topology_files_error_message(key),
                exception=InputFileError,
            )

    def _check_moldescriptor_needs_restart(self) -> None:
        """
        Raises
        ------
        InputFileError
            If the moldescriptor_file key is used without a restart file.
        """

        if self.moldescriptor_file is not None and self.restart_file is None:
            self.logger.error(
                (
                    "The moldescriptor_file key can only be "
                    "used in a meaningful way if a restart file is given."
                ),
                exception=InputFileError,
            )

    def _infer_required_topology_files(self) -> bool:
        """
        Infers missing topology files (the restart and moldescriptor file)
        that a key of the input file requires.
        """
        if self.restart_file is None:
            restart_file = self._infer_restart_file()

            if restart_file is not None:
                self._set_inferred_key(
                    self.restart_file_key,
                    restart_file,
                    "str",
                )
                self.logger.info(
                    f"Inferred restart_file='{restart_file}' from traj_files."
                )

        if self.moldescriptor_file is None:
            moldescriptor_file = self._infer_moldescriptor_file()

            if moldescriptor_file is not None:
                self._set_inferred_key(
                    self.moldescriptor_file_key,
                    moldescriptor_file,
                    "str",
                )
                self.logger.info(
                    (
                        "Inferred moldescriptor_file="
                        f"'{moldescriptor_file}' from traj_files."
                    )
                )

        return (
            self.restart_file is not None and
            self.moldescriptor_file is not None
        )

    def _infer_restart_file(self) -> str | None:
        """
        Infers the restart file from the first trajectory filename.
        """
        restart_file = self._restart_file_candidate()

        if restart_file is None:
            return None

        if Path(self.resolve_path(restart_file)).is_file():
            return restart_file

        return None

    def _infer_moldescriptor_file(self) -> str | None:
        """
        Infers the moldescriptor file from the first trajectory directory.
        """
        moldescriptor_file = self._moldescriptor_file_candidate()

        if moldescriptor_file is None:
            return None

        if Path(self.resolve_path(moldescriptor_file)).is_file():
            return moldescriptor_file

        return None

    def _restart_file_candidate(self) -> str | None:
        """
        Gets the restart file candidate inferred from the first trajectory.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return None

        return str(Path(traj_file).with_suffix(".rst"))

    def _moldescriptor_file_candidate(self) -> str | None:
        """
        Gets the moldescriptor file candidate inferred from the first trajectory.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return None

        return str(Path(traj_file).parent / self._moldescriptor_file_default)

    def _first_traj_file(self) -> str | None:
        """
        Gets the first trajectory file as written in the input file.

        The candidates derived from it are stored back into the input
        dictionary, so they must stay relative to the input file like
        every other filename there; the file properties resolve them.
        """
        traj_files = _parse_files(self.dictionary, self.traj_files_key)

        if not traj_files:
            return None

        return traj_files[0]

    def _missing_topology_files_error_message(self, key: str) -> str:
        """
        Gets the error message for missing topology files.

        Parameters
        ----------
        key : str
            The key (or keys) that require the topology files.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return (
                f"The {key} key requires both a restart file "
                "and a moldescriptor file. Could not infer missing files "
                "because no trajectory files were available."
            )

        return (
            f"The {key} key requires both a restart file "
            "and a moldescriptor file. Could not infer missing files "
            f"from trajectory file '{traj_file}'. "
            f"Tried restart_file='{self._restart_file_candidate()}' "
            "and moldescriptor_file="
            f"'{self._moldescriptor_file_candidate()}'."
        )

    def _set_inferred_key(self, key: str, value, value_type: str):
        """
        Sets a value in the input dictionary after it was inferred.
        """
        self.dictionary.dict[key] = (value, value_type, "inferred")
