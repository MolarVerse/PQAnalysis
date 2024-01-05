from ...io import PQAnalysisInputFileReader as Reader


class RDFInputFileReader(Reader):

    required_keys = [
        Reader.traj_files_key,
        Reader.reference_selection_key,
        Reader.target_selection_key,
        Reader.out_file_key
    ]

    known_keys = required_keys + [
        Reader.r_max_key,
        Reader.r_min_key,
        Reader.delta_r_key,
        Reader.n_bins_key,
        Reader.use_full_atom_info_key,
        Reader.log_file_key
    ]

    def __init__(self, filename: str) -> None:
        self.filename = filename
        super().__init__(filename)

    def read(self):
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.known_keys)

        if self.use_full_atom_info is None:
            self.use_full_atom_info = False
