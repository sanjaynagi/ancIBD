"""ancIBD: Identify IBD segments between pairs of individuals in ancient DNA."""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    ### Read the version from the installed package metadata, so that
    ### pyproject.toml stays the single place the version is declared.
    __version__ = _version("ancIBD")
except PackageNotFoundError:  # Running from a source tree, not an install
    __version__ = "unknown"
