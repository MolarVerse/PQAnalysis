from beartype.typing import List

from PQAnalysis.io import BaseWriter, FileWritingMode
from PQAnalysis.atomicSystem import AtomicSystem


class NEPWriter(BaseWriter):
    """
    A class to write NEP training and testing files.
    """

    def __init__(self,
                 filename: str,
                 mode: FileWritingMode | str = "w",
                 ) -> None:
        """
        Parameters
        ----------
        filename : str
            _description_
        mode : FileWritingMode | str, optional
            _description_, by default "w"
        """
        super().__init__(filename, mode)

    def write_from_atomic_system(self,
                                 atomic_system: AtomicSystem,
                                 use_coords,
                                 use_forces: bool = False,
                                 use_energy: bool = False,
                                 use_stress: bool = False,
                                 use_virial: bool = False,
                                 ) -> None:

        #  atomic_system: AtomicSystem = None,
        #  file_prefixes: List[str] | str = None,
        #  use_coords: bool = True,
        #  use_forces: bool = False,
        #  use_energy: bool = False,
        #  use_stress: bool = False,
        #  use_virial: bool = False,
        #  ) -> None:
