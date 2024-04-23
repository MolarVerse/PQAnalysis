"""
A module containing a class to read input files to setup the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` class.
"""

from __future__ import annotations

# local imports
from PQAnalysis.io import PQAnalysisInputFileReader as Reader


class DiffCalcInputFileReader(Reader):

    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc` class.
    """

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.traj_files_key,
        Reader.target_selection_key,
        Reader.window_size_key,
        Reader.gap_key,
        Reader.out_file_key,

    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        Reader.log_file_key,
        Reader.use_full_atom_info_key,
        Reader.restart_file_key,
        Reader.moldescriptor_file_key,
        Reader.n_start_key,
        Reader.n_stop_key,
        Reader.n_frame_key
    ]

    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            the filename of the input file
        """
        self.filename = filename
        super().__init__(filename)

    def read(self):
        """
        Reads the input file and parses it.
        It also sets the raw_input_file and the dictionary.
        It checks if all required keys are set and if all keys are known.

        Raises
        ------
        InputFileError
            if not all required keys are set in the input file
        InputFileWarning
            if unknown keys are set in the input file
        InputFileError
            if the no_intra_molecular key is used without a restart file
        InputFileError
            if the no_intra_molecular key is used without a moldescriptor file
        InputFileError
            if the moldescriptor_file key is used without a restart file
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)

        if self.no_intra_molecular is not None and (self.restart_file is None or self.moldescriptor_file is None):
            raise Reader.InputFileError(
                "The no_intra_molecular key can only be used if both a restart file and a moldescriptor file are given."
            )

        if self.moldescriptor_file is not None and self.restart_file is None:
            raise Reader.InputFileError(
                "The moldescriptor_file key can only be used in a meaningful way if a restart file is given."
            )


input_keys_documentation = f"""

For the diffcalc input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.traj_files_key}
        - The trajectory files to read. This can be a single file or a list of files.
    * - {Reader.target_selection_key}
        - The selection string to select the target atoms. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.reference atom.
    * - {Reader.window_size_key}
        - The window size of the running average.
    * - {Reader.gap_key}
        - The gap of the running average.
    * - {Reader.out_file_key}
        - The output file to write the diffcalc data to.


.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.log_file_key}
        - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
        - Whether to use full atom information for the selections.
    * - {Reader.restart_file_key}
        - The restart file to read the topology from.
    * - {Reader.moldescriptor_file_key}
        - The moldescriptor file to read the reference residues from.
    * - {Reader.n_start_key}    
        - The start frame of the diffcalc.
    * - {Reader.n_stop_key}
        - The stop frame of the diffcalc.
    * - {Reader.n_frame_key}
        - The number of frames to calculate the diffcalc.
        
Note
----
Optional keys does not mean that they are optional for the analysis. They are optional in the input file, but they might be required for the analysis. This means that if an optional keyword is specified other keywords might be required. For example:

- If the :code:`{Reader.moldescriptor_file_key}` key is specified, the :code:`{Reader.restart_file_key}` key is required in order use the reference residues in any meaningful way.
"""

DiffCalcInputFileReader.__doc__ += input_keys_documentation