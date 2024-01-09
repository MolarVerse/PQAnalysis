"""
A module containing a class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""

from ...io import PQAnalysisInputFileReader as Reader


class RDFInputFileReader(Reader):
    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
    """

    required_keys = [
        Reader.traj_files_key,
        Reader.reference_selection_key,
        Reader.target_selection_key,
        Reader.out_file_key,
        Reader.restart_file_key,
        Reader.moldescriptor_file_key,
    ]
    """List[str]: The required keys of the input file
    
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.traj_files_key`: The filenames of the trajectory files    
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.reference_selection_key`: The selection of the reference atoms. See also: :py:class:`~PQAnalysis.topology.selection.Selection`
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.target_selection_key`: The selection of the target atoms. See also: :py:class:`~PQAnalysis.topology.selection.Selection`
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.out_file_key`: The filename of the output file
    """

    optional_keys = required_keys + [
        Reader.r_max_key,
        Reader.r_min_key,
        Reader.delta_r_key,
        Reader.n_bins_key,
        Reader.use_full_atom_info_key,
        Reader.no_intra_molecular_key,
        Reader.log_file_key
    ]
    """List[str]: The optional keys of the input file
    
    | - :py:const:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.r_max_key`: The maximum radius of the RDF analysis in Angstrom
    | - :py:const:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.r_min_key`: The minimum radius of the RDF analysis in Angstrom
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.delta_r_key`: The width of the bins of the RDF analysis in Angstrom
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.n_bins_key`: The number of bins of the RDF analysis
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.use_full_atom_info_key`: If True, the full atom information is used for the selection of the reference and target atoms. If False only the element types without the atom names are used. This setting is only relevant if the selection is given with Atom objects. See also: :py:class:`~PQAnalysis.topology.selection.Selection`
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.no_intra_molecular`: If True, the intra molecular pairs are not considered in the RDF analysis. If False, the intra molecular pairs are considered. This setting is only relevant if the topology of the system contains residue information. See also: :py:class:`~PQAnalysis.topology.selection.Selection`, :py:class:`~PQAnalysis.topology.topology.Topology`
    | - :py:attr:`~PQAnalysis.io.inputFileReader.PQAnalysis.PQAnalysis_inputFileReader.PQAnalysisInputFileReader.log_file_key`: The filename of the log file
    """

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
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)
