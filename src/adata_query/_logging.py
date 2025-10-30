# -- import packages: ---------------------------------------------------------
import logging
import pathlib
import sys

# -- constants: ---------------------------------------------------------------
NAME = "adata_query"

# -- configure logger: --------------------------------------------------------
def configure_logging(
    name: str = NAME,
    log_file: str = f"{NAME}.log",
    log_dir: str = ".log_cache",
) -> logging.Logger:
    """Configure logging for the given name.

    Args:
        name (str): The name of the logger.
        log_file (str): The name of the log file.
        log_dir (str): The directory to save the log file.

    Returns:
        logger (logging.Logger)
    """
    logger = logging.getLogger(name)
    logger.setLevel(
        logging.DEBUG
    )  # Keep at DEBUG to allow file handler to capture all levels

    # Prevent adding handlers multiple times (important in notebooks!)
    if logger.hasHandlers():
        return logger

    # Create .log_cache directory if it doesn't exist
    log_dir = pathlib.Path(log_dir)
    log_dir.mkdir(exist_ok=True)

    # Update log_file path to include the directory
    log_file_path = log_dir / log_file

    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # File handler (DEBUG+)
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    stream_formatter = logging.Formatter(f"{NAME} [%(levelname)s]: %(message)s")

    # Stream handler (INFO+)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(stream_formatter)

    # Add filter to stream handler to only show INFO and above
    class InfoFilter(logging.Filter):
        def filter(self, record) -> bool:
            return record.levelno >= logging.INFO

    ch.addFilter(InfoFilter())

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
