"""Version helpers for source checkouts and installed builds."""

try:
    from ._version import __version__, __version_tuple__
except ModuleNotFoundError:
    __version__ = "0+unknown"
    __version_tuple__ = (0, "unknown")
