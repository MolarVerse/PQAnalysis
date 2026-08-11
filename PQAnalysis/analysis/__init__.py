"""
This is a collection of analysis subpackages.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .momentum import Momentum, check_momentum
    from .msd import MSD, MSDDiffusionFit, MSDInputFileReader, msd
    from .output import (
        AnalysisColumn,
        AnalysisOutputError,
        AnalysisOutputFormat,
        AnalysisSchema,
        AnalysisTable,
        convert_analysis_output,
        read_analysis_table,
    )
    from .rdf import RDF, RDFDataWriter, RDFInputFileReader, RDFLogWriter, rdf
    from .spectrum_broadening import broaden, build_spectrum
    from .vacf import VACF, VACFInputFileReader, vacf, vacf_spectrum
    from .vibrational import (
        VibrationalAnalysisInputFileReader,
        VibrationalAnalysisResult,
        vibrations,
    )

_EXPORTS = {
    "Momentum": ".momentum",
    "check_momentum": ".momentum",
    "AnalysisColumn": ".output",
    "AnalysisOutputError": ".output",
    "AnalysisOutputFormat": ".output",
    "AnalysisSchema": ".output",
    "AnalysisTable": ".output",
    "convert_analysis_output": ".output",
    "read_analysis_table": ".output",
    "MSD": ".msd",
    "MSDDiffusionFit": ".msd",
    "MSDInputFileReader": ".msd",
    "msd": ".msd",
    "RDF": ".rdf",
    "RDFInputFileReader": ".rdf",
    "RDFDataWriter": ".rdf",
    "RDFLogWriter": ".rdf",
    "rdf": ".rdf",
    "broaden": ".spectrum_broadening",
    "build_spectrum": ".spectrum_broadening",
    "VACF": ".vacf",
    "VACFInputFileReader": ".vacf",
    "vacf": ".vacf",
    "vacf_spectrum": ".vacf",
    "VibrationalAnalysisInputFileReader": ".vibrational",
    "VibrationalAnalysisResult": ".vibrational",
    "vibrations": ".vibrational",
}

__all__ = list(_EXPORTS)



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
