"""Progress-bar selection for terminal and notebook environments."""

import sys



def tqdm(*args, **kwargs):
    """Create a terminal bar unless the process is running in a notebook."""
    if "ipykernel" in sys.modules:
        from tqdm.auto import tqdm as implementation  # pylint: disable=import-outside-toplevel
    else:
        from tqdm import tqdm as implementation  # pylint: disable=import-outside-toplevel

    return implementation(*args, **kwargs)
