# -- configure logging: -------------------------------------------------------
from ._logging import configure_logging

logger = configure_logging()

# -- fetch version: -----------------------------------------------------------
# -- functional imports: ------------------------------------------------------
from . import _core
from .__version__ import __version__
from ._core import fetch, format_data, locate

# -- export: ------------------------------------------------------------------
__all__ = ["format_data", "fetch", "locate", "_core", "__version__", "logger"]
