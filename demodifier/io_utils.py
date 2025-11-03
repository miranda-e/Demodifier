# io_utils.py
# ---------------------------------------------------------
# Handles input file selection and reading for The Demodifier tool.
# Supports:
#   - Standard CSVs
#   - Mascot CSVs with preamble
#   - MaxQuant TSVs saved as .txt
# Includes a GUI file dialog for convenience.
# ---------------------------------------------------------

import os
import csv

def ask_for_csv_file():
    """
    Opens a file dialog to let the user select a CSV or TXT file
    (in Mascot or MaxQuant format). Uses a lazy import so it
    doesn't crash on headless systems.
    Returns
    -------
    str
        Selected file path, or an empty string if cancelled or unavailable.
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        # Headless environment or Tkinter not installed
        print("GUI file dialog not available (no display or tkinter missing).")
        return ""

    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_path = filedialog.askopenfilename(
        title="Select a CSV or TXT file",
        filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt")],
    )
    try:
        root.destroy()
    except Exception:
        pass
    return file_path or ""


def _iter_csv_with_preamble(file_obj):
    """
    Helper: yields dict rows from a CSV file that may have several rows of preamble
    of varying length followed by a table. The first header name for the table is
    always "prot_hit_num".

    Yields
    ------
    dict
        Each data row as an ordered dictionary.
    """
    header = None

    # Scan the file until the header row (starts with 'prot_hit_num') is found
    while True:
        pos = file_obj.tell()
        line = file_obj.readline()
        if not line:
            break
        cells = next(csv.reader([line]))
        if cells and cells[0].strip() == "prot_hit_num":
            header = cells
            break

    if header is None:
        # No preamble detected; reset and read normally
        file_obj.seek(0)
        reader = csv.DictReader(file_obj)
        yield from reader
        return
    # Yield the table starting from the detected header
    def line_iter():
        yield ",".join(header) + "\n"
        yield from file_obj

    reader = csv.DictReader(line_iter())
    yield from reader


def read_input_rows(path):
    """
    Unified row loader supporting:
      1) Standard CSV with header in the first row.
      2) CSV with several preamble rows; the real header row starts with 'prot_hit_num'.
      3) TSV saved with a .txt extension; header in the first row.

    Returns
    -------
    iterator of dict
        Iterator that yields each row as a dictionary.
    """
    _, ext = os.path.splitext(path)
    ext = ext.lower()

    if ext == ".txt":
        # MaxQuant format: tab-separated text files (.txt)
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            yield from reader
        return
    # Standard or Mascot-style CSV files
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        buf = fh.read(4096)
        fh.seek(0)
        if "prot_hit_num" in buf:
            yield from _iter_csv_with_preamble(fh)
        else:
            reader = csv.DictReader(fh)
            yield from reader
