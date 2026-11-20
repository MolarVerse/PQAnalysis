"""
Formatting helpers for self-describing analysis output files.
"""

import re
from collections.abc import Sequence

_HEADER_UNIT = re.compile(r"(?P<field>\S+)\s+\[(?P<unit>[^\]\s]+)\]")



def format_output_header(
    title: str,
    columns: Sequence[tuple[str, str, str]],
) -> str:
    """
    Build a compact metadata block describing a columnar output file.

    Parameters
    ----------
    title : str
        Human-readable name of the output.
    columns : Sequence[tuple[str, str, str]]
        ASCII field identifier, scientific Unicode symbol and unit for each
        column, in output order. Values must not contain whitespace.

    Returns
    -------
    str
        A title followed by machine-readable ``FIELDS``, human-readable
        ``SYMBOLS`` and ``UNITS`` comment lines. Standard numerical readers
        treat every line as a comment.
    """
    fields = " ".join(field for field, _, _ in columns)
    symbols = " ".join(symbol for _, symbol, _ in columns)
    units = " ".join(unit for _, _, unit in columns)

    return "\n".join(
        (
            f"# PQAnalysis: {title}",
            f"# FIELDS {fields}",
            f"# SYMBOLS {symbols}",
            f"# UNITS {units}",
        )
    )



def header_name(field: str, unit: str, with_unit: bool) -> str:
    """Return a CSV/TSV header entry, ``FIELD [unit]`` if requested."""
    if not with_unit or unit == "1":
        return field
    return f"{field} [{unit}]"



def split_header_units(
    header: Sequence[str],
) -> tuple[tuple[str, ...], tuple[str | None, ...]]:
    """
    Split CSV/TSV header entries into fields and optional units.

    >>> split_header_units(["E(TOT) [kcal/mol]", "N"])
    (('E(TOT)', 'N'), ('kcal/mol', None))
    """
    fields, units = [], []
    for value in header:
        match = _HEADER_UNIT.fullmatch(value.strip())
        fields.append(match["field"] if match else value.strip())
        units.append(match["unit"] if match else None)
    return tuple(fields), tuple(units)
