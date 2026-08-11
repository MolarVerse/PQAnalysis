"""
A package containing classes and functions with common use.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .common import __header__, print_header
    from .decorators import (
        count_decorator,
        instance_function_count_decorator,
        timeit_in_class,
    )
    from .units import calculate_simulation_time

_EXPORTS = {
    "print_header": ".common",
    "__header__": ".common",
    "count_decorator": ".decorators",
    "instance_function_count_decorator": ".decorators",
    "timeit_in_class": ".decorators",
    "calculate_simulation_time": ".units",
}

__all__ = [name for name in _EXPORTS if not name.startswith("_")]



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
