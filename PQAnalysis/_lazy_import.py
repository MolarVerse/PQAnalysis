"""Utilities for preserving package exports without eager imports."""

from importlib import import_module
from types import ModuleType
from typing import Any



def resolve_export(
    package_name: str,
    namespace: dict[str, Any],
    exports: dict[str, str],
    name: str,
) -> Any:
    """Import and cache one public package attribute on first access."""
    try:
        module_name = exports[name]
    except KeyError as exception:
        raise AttributeError(
            f"module {package_name!r} has no attribute {name!r}"
        ) from exception

    value = getattr(import_module(module_name, package_name), name)

    child_name = module_name.lstrip(".").partition(".")[0]
    if (
        child_name != name and child_name in exports and
        isinstance(namespace.get(child_name), ModuleType)
    ):
        namespace.pop(child_name)

    namespace[name] = value
    return value



def public_dir(
    namespace: dict[str, Any],
    exports: dict[str, str],
) -> list[str]:
    """Return loaded and deferred package attributes for introspection."""
    return sorted(set(namespace) | set(exports))
