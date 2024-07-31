"""
This module provides API functions for bulk modulus analysis.
"""
import numpy as np
from beartype.typing import List, Tuple

from PQAnalysis.io.energy_file_reader import EnergyFileReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking

from PQAnalysis.types import Np1DNumberArray, Np2DNumberArray
from .bulk_modulus import BulkModulus
from .bulk_modulus_input_file_reader import BulkModulusInputFileReader
from .bulk_modulus_output_file_writer import BulkModulusDataWriter
from .bulk_modulus_output_file_writer import BulkModulusLogWriter


@runtime_type_checking
def bulk_modulus(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the bulk modulus using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit
    :py:mod:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus_input_file_reader`.
    For more information on the exact calculation of the bulk modulus
    :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit
        :py:mod:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """
    md_format = MDEngineFormat(md_format)

    input_reader = BulkModulusInputFileReader(input_file)
    input_reader.read()

    if input_reader.energy_perturbation_files is None:
        raise ValueError(
            "The energy perturbation files are required."
        )

    if md_format != MDEngineFormat.PQ:
        use_info_file = True
        info_file = input_reader.info_file

    else:
        info_file = input_reader.info_file
        use_info_file = False

    if (len(input_reader.energy_perturbation_files) == 0):
        raise ValueError(
            "The energy perturbation files are required."
        )
    if (len(input_reader.volumes_perturbation) == 0):
        raise ValueError(
            "The volumes perturbation are required."
        )

    if (len(input_reader.energy_perturbation_files) != len(input_reader.volumes_perturbation)):
        raise ValueError(
            "The number of energy perturbation files and volumes perturbation must be the same."
        )

    if (len(input_reader.energy_perturbation_files) == 1):
        raise ValueError(
            "At least two energy perturbation files are required."
        )

    if (len(input_reader.volumes_perturbation) == 1):
        raise ValueError(
            "At least two volumes perturbation are required."

        )

    pressures_perturbation = np.empty(
        (len(input_reader.volumes_perturbation), 2))

    for i, file in enumerate(input_reader.energy_perturbation_files):

        energy_reader = EnergyFileReader(
            filename=file,
            info_filename=info_file,
            use_info_file=use_info_file,
            engine_format=md_format)
        energy_data = energy_reader.read()
        pressure_units = energy_data.units["PRESSURE"]
        pressure_data = energy_data.data[energy_data.info["PRESSURE"]]
        pressures_perturbation[i, 0] = np.average(pressure_data)
        pressures_perturbation[i, 1] = np.std(pressure_data)

    if input_reader.volumes_perturbation is None:
        raise ValueError(
            "The volumes perturbation must be provided."
        )

    if input_reader.volume_equilibrium is None:
        raise ValueError(
            "The equilibrium volume must be provided."
        )

    if input_reader.mode is None:
        mode = "simple"
    else:
        mode = input_reader.mode

    _bulk_modulus = BulkModulus(
        pressures_perturbation=pressures_perturbation,
        volumes_perturbation=np.array(input_reader.volumes_perturbation),
        volume_equilibrium=input_reader.volume_equilibrium,
        mode=mode
    )

    data_writer = BulkModulusDataWriter(
        filename=input_reader.out_file_key
    )

    log_writer = BulkModulusLogWriter(filename=input_reader.log_file)

    log_writer.write_before_run(_bulk_modulus)

    bulk_modulus_data = _bulk_modulus.run()

    data_writer.write(Tuple(bulk_modulus_data))
    log_writer.write_after_run(_bulk_modulus, pressure_units)
