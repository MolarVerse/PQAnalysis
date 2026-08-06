"""
Formatting helpers for self-describing analysis output files.
"""

from collections.abc import Sequence



def format_output_header(
    title: str,
    columns: Sequence[tuple[str, str]],
) -> str:
    """
    Build a compact metadata block describing a columnar output file.

    Parameters
    ----------
    title : str
        Human-readable name of the output.
    columns : Sequence[tuple[str, str]]
        Field symbol and unit for each column, in output order. Symbols and
        units must not contain whitespace.

    Returns
    -------
    str
        A title followed by machine-readable ``FIELDS`` and ``UNITS`` comment
        lines. Standard numerical readers treat every line as a comment.
    """
    fields = " ".join(field for field, _ in columns)
    units = " ".join(unit for _, unit in columns)

    return "\n".join(
        (
            f"# PQAnalysis: {title}",
            f"# FIELDS {fields}",
            f"# UNITS {units}",
        )
    )
