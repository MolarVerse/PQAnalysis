"""
This module provides API functions for the velocity auto-correlation
function (VACF) analysis.
"""

import numpy as np

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.type_checking import runtime_type_checking

from .vacf import VACF
from .spectrum import vacf_spectrum
from .vacf_input_file_reader import VACFInputFileReader
from .vacf_output_file_writer import (
    VACFDataWriter,
    VACFLogWriter,
    VACFSpectrumDataWriter,
    VACFWindowedDataWriter,
)
from .exceptions import VACFError

#: The default correlation window size in frames.
WINDOW_SIZE_DEFAULT = 1000



@runtime_type_checking
def read_static_charges(
    filename: str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
) -> Np1DNumberArray:
    """
    Reads a legacy static charge file.

    The file consists of a header line with the number of atoms, one
    comment line and one "name charge" line per atom (legacy
    Fluxfreqcalc chrg_file format). For the
    :py:class:`~PQAnalysis.traj.formats.MDEngineFormat.QMCFC` format
    the first entry is the dummy 'X' atom, which is stripped like in
    the velocity trajectory.

    Parameters
    ----------
    filename : str
        The static charge file to read.
    md_format : MDEngineFormat | str, optional
        the format of the underlying trajectory. Default is "PQ".

    Returns
    -------
    Np1DNumberArray
        The static atomic partial charges.

    Raises
    ------
    VACFError
        If the file header or a charge line cannot be parsed, if the
        file provides fewer charges than the header announces or if
        the first atom of a QMCFC charge file is not the dummy 'X'
        atom.
    """
    md_format = MDEngineFormat(md_format)

    with open(filename, "r", encoding="utf-8") as file:
        lines = file.readlines()

    try:
        n_entries = int(lines[0].split()[0])
    except (IndexError, ValueError) as exception:
        raise VACFError(
            f"Could not read the number of atoms from the header of "
            f"the charge file '{filename}'."
        ) from exception

    # lines[1] is a comment line and is skipped like in the legacy tool
    entry_lines = lines[2:2 + n_entries]

    if len(entry_lines) != n_entries:
        raise VACFError(
            f"The charge file '{filename}' provides only "
            f"{len(entry_lines)} charge line(s), but the header "
            f"announces {n_entries} atom(s)."
        )

    names = []
    charges = []

    for line in entry_lines:
        splitted_line = line.split()

        try:
            names.append(splitted_line[0])
            charges.append(float(splitted_line[1]))
        except (IndexError, ValueError) as exception:
            raise VACFError(
                f"Could not parse the charge line '{line.strip()}' of "
                f"the charge file '{filename}'. Each line has to be "
                "of the form 'name charge'."
            ) from exception

    if md_format == MDEngineFormat.QMCFC:
        if not names or names[0].upper() != "X":
            raise VACFError(
                f"The first atom of the QMCFC charge file '{filename}' "
                "is not the dummy 'X' atom."
            )

        charges = charges[1:]

    return np.array(charges, dtype=np.float64)



@runtime_type_checking
def vacf(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the velocity auto-correlation function (VACF) using a
    given input file.

    This is just a wrapper function combining the underlying classes
    and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.vacf.vacf_input_file_reader`.
    For more information on the exact calculation of
    the VACF please visit :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF`
    and for the spectrum
    :py:func:`~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit
        :py:mod:`~PQAnalysis.analysis.vacf.vacf_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """

    md_format = MDEngineFormat(md_format)

    input_reader = VACFInputFileReader(input_file)
    input_reader.read()

    traj_reader = TrajectoryReader(
        input_reader.traj_files,
        md_format=md_format
    )

    charges = None
    charge_traj = None

    if input_reader.charge_file is not None:
        charges = read_static_charges(
            input_reader.charge_file,
            md_format=md_format,
        )

    if input_reader.charge_files is not None:
        charge_traj = TrajectoryReader(
            input_reader.charge_files,
            traj_format=TrajectoryFormat.CHARGE,
            md_format=md_format,
        )

    window_size = input_reader.window

    if window_size is None:
        window_size = WINDOW_SIZE_DEFAULT

    _vacf = VACF(
        traj=traj_reader,
        window_size=window_size,
        time_step=input_reader.time_step,
        target_species=input_reader.target_selection,
        gap=input_reader.gap,
        charges=charges,
        charge_traj=charge_traj,
        method=input_reader.method,
        use_full_atom_info=input_reader.use_full_atom_info,
    )

    # all output writers are constructed before the run so that a
    # pre-existing output file aborts before the expensive analysis
    data_writer = VACFDataWriter(input_reader.out_file)
    log_writer = VACFLogWriter(input_reader.log_file)

    spectrum_writer = None
    windowed_writer = None

    if input_reader.spectrum_file is not None:
        spectrum_writer = VACFSpectrumDataWriter(input_reader.spectrum_file)

        if input_reader.windowed_out_file is not None:
            windowed_writer = VACFWindowedDataWriter(
                input_reader.windowed_out_file
            )

    log_writer.write_before_run(_vacf)

    time, correlation = _vacf.run()

    data_writer.write((time, correlation))

    if spectrum_writer is not None:
        spectrum_kwargs = {}

        if input_reader.ftsize is not None:
            spectrum_kwargs["ftsize"] = input_reader.ftsize
        if input_reader.window_function is not None:
            spectrum_kwargs["window_function"] = input_reader.window_function
        if input_reader.window_param is not None:
            spectrum_kwargs["window_param"] = input_reader.window_param
        if input_reader.window_start is not None:
            spectrum_kwargs["window_start"] = input_reader.window_start
        if input_reader.window_stop is not None:
            spectrum_kwargs["window_stop"] = input_reader.window_stop

        wavenumbers, amplitudes, windowed_correlation = vacf_spectrum(
            time,
            correlation,
            **spectrum_kwargs,
        )

        spectrum_writer.write((wavenumbers, amplitudes))

        if windowed_writer is not None:
            windowed_writer.write((time, windowed_correlation))

    log_writer.write_after_run(_vacf)
