# settings.py
# ---------------------------------------------------------
# Handles:
#   - Logging setup for the entire tool
#   - Command-line user prompts for threads / verbosity
# ---------------------------------------------------------

import os
import logging

# ---------------------------------------------------------
# Create a clean logger (no side effects on import)
# ---------------------------------------------------------
logger = logging.getLogger("demodifier")
# No handler attached yet; configured later via setup_logging()


# ---------------------------------------------------------
# CLI settings helper
# ---------------------------------------------------------
def ask_for_settings_cli(default_threads=None, default_verbose=False):
    """
    Ask the user in the terminal:
      How many processors?
      Verbose mode on?
    Returns (num_threads:int, verbose:bool).
    """

    # Default number of processors is system CPU count or 4
    if default_threads is None:
        try:
            default_threads = max(1, (os.cpu_count() or 4))
        except Exception:
            default_threads = 4

    # Prompt for processors
    while True:
        try:
            print("How many processors?")
            raw = input().strip()
            num_threads = int(raw)
            if num_threads < 1:
                raise ValueError
            break
        except Exception:
            print("Please enter a whole number ≥ 1.")

    # Prompt for verbose
    while True:
        print("Verbose mode on?")
        raw = input().strip().lower()
        if raw in {"y", "yes"}:
            verbose = True
            break
        if raw in {"n", "no"}:
            verbose = False
            break
        if raw == "":
            verbose = default_verbose
            break
        print("Please answer yes or no.")

    return num_threads, verbose


# ---------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------
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