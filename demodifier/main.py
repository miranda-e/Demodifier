# main.py
# ---------------------------------------------------------
# Orchestrates the Demodifier pipeline:
#   1) Load input rows (CSV/TSV, including Mascot preamble handling).
#   2) Generate permutations + query Unipept for LCAs in parallel.
#   3) Delegate all writing to writers.write_results (CSV + JSON).
#
# CLI behavior:
#   • If you pass an input path: no GUI dialog is used for file selection.
#   • If you don't pass a path: we try to open a file GUI dialog.
#   • If the dialog isn't available (headless) or is canceled: we print guidance and exit.
#   • --threads / --verbose [yes|no] override prompts. Bare --verbose = yes.
#   • If --threads is omitted, you will be PROMPTED for a value.
#   • If --verbose is omitted, you will be PROMPTED (bare --verbose = yes).
#   • When using bare --verbose, place it after the input file (positional must come first).
# ---------------------------------------------------------

import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from .analysis import (
    extract_deamidation_count,
    generate_deamidation_permutations,
    generate_reamidation_permutations,
    add_pyro_glu_permutations,
)
from .unipept_api import process_peptides, get_lcas_for_permutations, make_session
from .io_utils import ask_for_csv_file, read_input_rows
from .writers import write_results
from .settings import logger, setup_logging

# Per-thread HTTP session (each worker keeps its own keep-alive pool)
_thread_local = threading.local()

# Each worker thread gets its own HTTP session for stability (avoids cross-thread sharing).
def _get_thread_session():
    if not hasattr(_thread_local, "session"):
        _thread_local.session = make_session()
    return _thread_local.session

def process_row(row):
    """
    Process one input row:
      - Parse fields
      - Generate MISPs
      - Query Unipept for LCA(s)
    Returns (row, final_permutations, row_lcas, input_pep_lca) or None if no peptide.
    """
    session = _get_thread_session()

    peptide = (row.get('Sequence') or row.get('pep_seq') or '').strip()
    modifications = (row.get('Modifications') or row.get('pep_var_mod') or '').strip()

    logger.debug(f"Processing row: Peptide={peptide!r}, Mods={modifications!r}")
    if not peptide:
        return None

    max_substitutions = extract_deamidation_count(modifications)
    peptide_options = generate_deamidation_permutations(peptide, max_substitutions, modifications)

    final_permutations = []
    for perm, _, modified_positions in peptide_options:
        final_permutations.extend(
            generate_reamidation_permutations(perm, modified_positions, modifications)
        )

    # Add optional pyro-Glu-derived permutations; then deduplicate
    final_permutations = list(set(add_pyro_glu_permutations(final_permutations, modifications)))

    input_pep_lca = process_peptides([peptide], session=session)[0]
    row_lcas = get_lcas_for_permutations(final_permutations, session)

    return row, final_permutations, row_lcas, input_pep_lca


def make_output_paths(input_csv: str):
    """
    Build output file paths alongside the input.
    """
    base = os.path.splitext(input_csv)[0]
    return (
        f"{base}_results.csv",        # summary CSV
        f"{base}_output.json",        # JSON summary
        f"{base}_permutations.csv",   # per-variant CSV
    )


def run_pipeline(input_csv: str, num_threads: int):
    """
    Process all rows in a thread pool and collect results in memory.
    """
    results = []
    rows_iter = read_input_rows(input_csv)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_row, row) for row in rows_iter]
        for future in as_completed(futures):
            result = future.result()
            if result is None:
                continue
            results.append(result)

    return results


def main(input_csv=None, num_threads=None, verbose=None):
    """
    CLI entry.
    - If no path: try a GUI file dialog.
    - If the dialog is unavailable or canceled: print guidance and exit.
    - If --threads is omitted: prompt “How many processors?”
    - If --verbose is omitted: prompt “Verbose mode on? (yes/no)” (bare --verbose = yes).
    - Configure logging, run pipeline, write results.
    """
    # 1) Resolve input path
    if input_csv is None:
        input_csv = ask_for_csv_file()

    if not input_csv:
        # Dialog canceled OR not available (e.g., headless environment)
        print(
            "No file selected and GUI picker unavailable.\n"
            "Please run again and provide an input file path, for example:\n"
            "    python -m demodifier.main\n"
            "Optionally add flags, e.g.:\n"
            "    demodifier your_input.csv --threads 8 --verbose yes"
        )
        return

    # 2) Resolve execution settings (PROMPT when flags not provided)
    if num_threads is None:
        while True:
            try:
                raw = input("How many processors?\n").strip()
                num_threads = int(raw)
                if num_threads < 1:
                    raise ValueError
                break
            except Exception:
                print("Please enter a whole number ≥ 1.")

    if verbose is None:
        # Ask only if not provided via --verbose
        while True:
            raw = input("Verbose mode on? (yes/no)\n").strip().lower()
            if raw in {"y", "yes"}:
                verbose = True
                break
            if raw in {"n", "no"}:
                verbose = False
                break
            print("Please answer yes or no.")

    # 3) Configure logging
    setup_logging(verbose)

    # 4) Run pipeline and write outputs
    start = time.time()
    output_csv, output_json, output_variant_lca_csv = make_output_paths(input_csv)
    logger.debug(f"Starting Demodifier with {num_threads} threads on {input_csv!r}...")

    results = run_pipeline(input_csv, num_threads)
    write_results(results, output_csv, output_json, output_variant_lca_csv)

    elapsed = time.time() - start
    logger.debug(f"Completed in {elapsed:.2f} seconds with {num_threads} threads.")
    print("\nDemodifier completed")

# --- Command-line interface ---
# Standard Python entry point — runs when script is executed directly.
# Adds optional flags that allow users to skip prompts when flags are provided.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process peptides from a CSV/TSV file.")
    parser.add_argument("input_csv", nargs="?", default=None, help="Path to input CSV/TSV file")
    parser.add_argument("--threads", type=int, help="Number of worker threads (skip prompt)")
    # Accepts optional value (yes/no). Bare --verbose => True.
    # Any value not in {1,true,t,y,yes,on} is treated as False; order matters if you use the bare flag.
    parser.add_argument(
        "--verbose",
        nargs="?",              # optional value allowed
        const=True,             # bare --verbose => True
        default=None,           # omitted => prompt later
        type=lambda x: str(x).lower() in ("1", "true", "t", "y", "yes", "on"),
        help=(
            "Enable verbose logging (yes/no). Bare --verbose = yes. "
            "Note: place --verbose after the input file if using the bare flag, e.g. "
            "'demodifier input.csv --threads 4 --verbose'."
        ),
    )
    args = parser.parse_args()

    # Warn if the arguments look out of order (e.g., '--verbose example.csv ...')
    argv = os.sys.argv[1:]
    if args.input_csv is None and "--verbose" in argv:
        idx = argv.index("--verbose")
        if idx + 1 < len(argv):
            next_tok = argv[idx + 1]
            if not next_tok.startswith("-"):
                print("\nIt looks like your arguments might be in the wrong order.")
                print("file path must be first argument.")
                print("Example of correct usage:")
                print("    demodifier your_input.csv --verbose yes --threads 4")
                print("Or:")
                print("    demodifier your_input.csv --threads 4 --verbose\n")

    main(
        input_csv=args.input_csv,
        num_threads=args.threads if args.threads is not None else None,
        verbose=args.verbose,
    )

