class _OutputFileMixin:
    @property
    def start_file(self):
        return self.dictionary["start_file"][0]

    @property
    def rpmd_start_file(self):
        return self.dictionary["rpmd_start_file"][0]

    @property
    def restart_file(self):
        return self.dictionary["rst_file"][0]

    @property
    def trajectory_file(self):
        return self.dictionary["traj_file"][0]

    @property
    def velocity_file(self):
        return self.dictionary["vel_file"][0]

    @property
    def force_file(self):
        return self.dictionary["force_file"][0]

    @property
    def charge_file(self):
        return self.dictionary["charge_file"][0]

    @property
    def energy_file(self):
        return self.dictionary["energy_file"][0]

    @property
    def info_file(self):
        return self.dictionary["info_file"][0]

    @property
    def output_file(self):
        return self.dictionary["output_file"][0]

    @property
    def file_prefix(self):
        return self.dictionary["file_prefix"][0]

    @property
    def rpmd_restart_file(self):
        return self.dictionary["rpmd_rst_file"][0]

    @property
    def rpmd_trajectory_file(self):
        return self.dictionary["rpmd_traj_file"][0]

    @property
    def rpmd_velocity_file(self):
        return self.dictionary["rpmd_vel_file"][0]

    @property
    def rpmd_force_file(self):
        return self.dictionary["rpmd_force_file"][0]

    @property
    def rpmd_energy_file(self):
        return self.dictionary["rpmd_energy_file"][0]

    @property
    def rpmd_charge_file(self):
        return self.dictionary["rpmd_charge_file"][0]
