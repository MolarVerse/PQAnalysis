"""
A module containing functions to work with files, which are of general use.
"""

from pathlib import Path

import glob
import numpy as np

from beartype.typing import List



def find_files_with_prefix(file_prefixes: List[str] | str) -> List[str]:
    """
    Finds files with a given prefix.

    Parameters
    ----------
    file_prefixes : List[str] or str
        The prefixes of the files to find. Here with prefix we mean
        the part of the filename not only the name before the extension,
        but every matching file that starts with the given prefix.

    Returns
    -------
    List[str]
        The files with the given prefixes.
    """

    file_prefixes = list(np.atleast_1d(file_prefixes))
    file_prefixes = [file_prefix + "*" for file_prefix in file_prefixes]
    file_prefixes = [glob.glob(prefix + "*") for prefix in file_prefixes]

    files = [
        file for sublist in file_prefixes for file in sublist
        if Path(file).is_file()
    ]

    return files
