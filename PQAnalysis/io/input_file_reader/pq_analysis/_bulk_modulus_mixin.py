"""
A module containing a Mixin class for bulk modulus 
related keywords of a PQAnalysis input file.
"""


from ._parse import _parse_float_list, _parse_real


class _BulkModulusMixin:

    """
    A mixin class to read all file related
    keywords from the input dictionary.

    The following keywords are read:
        - volumes_perturbation
        - volumes_equilibrium
    """

    @property
    def volumes_perturbation(self) -> list[float] | None:
        """list[float] | None: The perturbation volumes of the simulation."""
        return _parse_float_list(self.dictionary, self.volumes_perturbation_key)

    @property
    def volume_equilibrium(self) -> float | None:
        """float | None: The equilibrium volume of the simulation."""
        return _parse_real(self.dictionary, self.volume_equilibrium_key)
