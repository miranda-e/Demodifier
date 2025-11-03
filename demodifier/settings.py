# settings.py
# ---------------------------------------------------------
# Handles:
#   - Logging setup for the entire tool
#   - Command-line user prompts for threads / verbosity

import os
import logging

# Create a clean logger for the Demodifier tool
logger = logging.getLogger("demodifier")

# Logging configuration
def setup_logging(verbose: bool):
    """
    Configure logging once, without side effects on import.
    Attaches a console handler only if none exists.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(level)

    # Avoid duplicate logs if setup_logging() is called multiple times
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Prevent messages from being duplicated through root logger
    logger.propagate = False
    return logger