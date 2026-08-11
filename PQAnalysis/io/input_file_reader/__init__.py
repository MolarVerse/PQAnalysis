"""
A package containing modules, classes and functions to 
parse and read input files of PQAnalysis itself and the md engines it supports.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .formats import InputFileFormat
    from .input_file_parser import InputDictionary, InputFileParser
    from .pq.pq_input_file_reader import PQInputFileReader
    from .pq_analysis.pqanalysis_input_file_reader import (
        PQAnalysisInputFileReader,
    )

_EXPORTS = {
    "InputFileParser": ".input_file_parser",
    "InputDictionary": ".input_file_parser",
    "PQInputFileReader": ".pq.pq_input_file_reader",
    "PQAnalysisInputFileReader": ".pq_analysis.pqanalysis_input_file_reader",
    "InputFileFormat": ".formats",
}

__all__ = list(_EXPORTS)



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
