# -- import packages: ---------------------------------------------------------
import importlib.metadata

# -- function: ----------------------------------------------------------------
def fetch_version():
    try:
        __version__ = importlib.metadata.version("adata-query")
    except importlib.metadata.PackageNotFoundError:
        # Package is not installed, so we can't determine the version
        __version__ = "unknown"
    return __version__

# -- fetch version: -----------------------------------------------------------
__version__ = fetch_version()

# -- export: ------------------------------------------------------------------
__all__ = ["__version__"]