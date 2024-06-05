"""
This module contains utility functions for string manipulation.
"""



def is_comment_line(line: str, comment_token="#", empty_line=False) -> bool:
    """
    Returns whether the line is a comment line.

    Parameters
    ----------
    line : str
        The line to be checked.
    comment_token : str, optional
        The token that indicates a comment line, by default "#".
    empty_line : bool, optional
        Whether an empty line should be considered a comment line, by default False.

    Returns
    -------
    bool
        Whether the line is a comment line.
    """

    if empty_line and not line.strip():
        return True

    return line.strip().startswith(comment_token)
