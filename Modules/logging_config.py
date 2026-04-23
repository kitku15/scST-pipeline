import logging
import time
from pathlib import Path


def configure_logging(logging_dir: Path | None = None, log_level=logging.INFO):
    """Configure logging for the package.

    Args:
        logging_dir: The directory where log files will be stored.
        log_level: The logging level to set. Defaults to logging.INFO.
    """
    if logging_dir is None:
        logging_dir = Path("analysis/logs")

    logger = (
        logging.getLogger()
    )  # Configure root logger so all loggers inherit handlers
    logger.setLevel(log_level)

    ch = logging.StreamHandler()
    fh = logging.FileHandler(
        logging_dir / f"recode_st-{time.strftime('%Y%m%d-%H%M%S')}.log"
    )

    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] %(message)s",
        datefmt="%Y/%m/%d %H:%M:%S",
    )

    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)
