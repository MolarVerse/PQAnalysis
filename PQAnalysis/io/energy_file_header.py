"""
Column names and units from the header of PQ energy files.

The header starts the ``.en`` file and ends with the column line
directly above the data::

    # PQ energy file v1.1.0
    # timestep = 0.5 fs
    #STEP[-] TEMPERATURE[K] PRESSURE[bar] ...

The column line has one ``NAME[unit]`` token per data column; comment
lines such as the timestep may sit between the first and the last line.
The step column ``STEP[-]`` is read as the older
``SIMULATION-TIME[step]``.
"""

import re

from .exceptions import EnergyFileReaderError

EnergyLabels = tuple[dict[str, int], dict[str, str]]

_HEADER_PATTERN = re.compile(r"#\s*PQ energy file(\s|$)")
_COLUMN_PATTERN = re.compile(r"(?P<name>[^\[\]\s]+)\[(?P<unit>[^\[\]\s]*)\]")
_COLUMN_ALIASES = {("STEP", "-"): ("SIMULATION-TIME", "step")}



def read_energy_header(filename: str) -> EnergyLabels | None:
    """
    Read the column names and units from a PQ energy-file header.

    Parameters
    ----------
    filename : str
        Energy file.

    Returns
    -------
    tuple[dict[str, int], dict[str, str]] | None
        Column index and unit per name, in the form of
        :py:meth:`~PQAnalysis.io.InfoFileReader.read`, or None if the
        first line is not ``# PQ energy file``.

    Raises
    ------
    EnergyFileReaderError
        If the column line is missing or malformed.
    """
    with open(filename, encoding="utf-8") as file:
        if _HEADER_PATTERN.match(file.readline().strip()) is None:
            return None
        column_line = _last_comment_line(file)

    if column_line is None or not column_line[1:].split():
        raise EnergyFileReaderError(
            f"The PQ energy file header of {filename} needs a column line "
            "'#NAME[unit] NAME[unit] ...' above the data."
        )

    info = {}
    units = {}
    for index, token in enumerate(column_line[1:].split()):
        match = _COLUMN_PATTERN.fullmatch(token)
        if match is None:
            raise EnergyFileReaderError(
                f"Invalid column '{token}' in the header of {filename}. "
                "Expected NAME[unit]."
            )
        name, unit = _COLUMN_ALIASES.get(
            (match["name"], match["unit"]), (match["name"], match["unit"])
        )
        if name in info:
            raise EnergyFileReaderError(
                f"Duplicate column '{match['name']}' in the header of "
                f"{filename}."
            )
        info[name] = index
        units[name] = unit

    return info, units



def _last_comment_line(file) -> str | None:
    last = None
    for line in file:
        text = line.strip()
        if text and not text.startswith("#"):
            break
        if text:
            last = text
    return last
