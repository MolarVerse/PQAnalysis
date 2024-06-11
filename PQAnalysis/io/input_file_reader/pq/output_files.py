"""
A module containing a Mixin class for output files of a PQ input file.
"""



class _OutputFileMixin:

    """
    Mixin class for output files of a PQ input file.
    """

    start_file_keys = ["start_file", "rpmd_start_file"]

    output_file_keys = [
        "restart_file",
        "traj_file",
        "vel_file",
        "force_file",
        "charge_file",
        "energy_file",
        "info_file",
        "output_file",
        "file_prefix",
        "rpmd_restart_file",
        "rpmd_traj_file",
        "rpmd_vel_file",
        "rpmd_force_file",
        "rpmd_energy_file",
        "rpmd_charge_file"
    ]

    @property
    def start_file(self) -> str:
        """str: The start file of the simulation."""
        return self.dictionary["start_file"][0]

    @property
    def is_start_file_defined(self) -> bool:
        """bool: Whether a start file is defined in the input file."""

        if "start_file" in self.dictionary.keys():
            return True

        return False

    @property
    def rpmd_start_file(self) -> str:
        """str: The rpmd start file of the simulation."""
        return self.dictionary["rpmd_start_file"][0]

    @property
    def is_rpmd_start_file_defined(self) -> bool:
        """bool: Whether a rpmd start file is defined in the input file."""
        if "rpmd_start_file" in self.dictionary.keys():
            return True

        return False

    @property
    def restart_file(self) -> str:
        """str: The restart file of the simulation."""
        return self.dictionary["restart_file"][0]

    @property
    def trajectory_file(self) -> str:
        """str: The trajectory file of the simulation."""
        return self.dictionary["traj_file"][0]

    @property
    def velocity_file(self) -> str:
        """str: The velocity file of the simulation."""
        return self.dictionary["vel_file"][0]

    @property
    def force_file(self) -> str:
        """str: The force file of the simulation."""
        return self.dictionary["force_file"][0]

    @property
    def charge_file(self) -> str:
        """str: The charge file of the simulation."""
        return self.dictionary["charge_file"][0]

    @property
    def energy_file(self) -> str:
        """str: The energy file of the simulation."""
        return self.dictionary["energy_file"][0]

    @property
    def info_file(self) -> str:
        """str: The info file of the simulation."""
        return self.dictionary["info_file"][0]

    @property
    def output_file(self) -> str:
        """str: The output file of the simulation."""
        return self.dictionary["output_file"][0]

    @property
    def file_prefix(self) -> str:
        """str: The file prefix of the simulation."""
        return self.dictionary["file_prefix"][0]

    @property
    def rpmd_restart_file(self) -> str:
        """str: The rpmd restart file of the simulation."""
        return self.dictionary["rpmd_restart_file"][0]

    @property
    def rpmd_trajectory_file(self) -> str:
        """str: The rpmd trajectory file of the simulation."""
        return self.dictionary["rpmd_traj_file"][0]

    @property
    def rpmd_velocity_file(self) -> str:
        """str: The rpmd velocity file of the simulation."""
        return self.dictionary["rpmd_vel_file"][0]

    @property
    def rpmd_force_file(self) -> str:
        """str: The rpmd force file of the simulation."""
        return self.dictionary["rpmd_force_file"][0]

    @property
    def rpmd_energy_file(self) -> str:
        """str: The rpmd energy file of the simulation."""
        return self.dictionary["rpmd_energy_file"][0]

    @property
    def rpmd_charge_file(self) -> str:
        """str: The rpmd charge file of the simulation."""
        return self.dictionary["rpmd_charge_file"][0]
