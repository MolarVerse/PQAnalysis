from beartype.typing import List, Tuple, Any

from ..exceptions import InputFileError


def parse_traj_files_key(data: Tuple[Any, str, str]) -> List[str]:
    type = data[1]

    if type == "str":
        return [data[0]]
    elif type == "glob" or type == "list(str)":
        return data[0]
    else:
        raise InputFileError(
            f"The traj_files key word has to be either a string, glob or a list of strings - actual it is parsed as a {type}")
