# main.py
# ---------------------------------------------------------
# Orchestrates the Demodifier pipeline:
#   1) Load input rows (CSV/TSV, including Mascot preamble handling).
#   2) Generate permutations + query Unipept for LCAs in parallel.
#   3) Delegate all writing to writers.write_results (CSV + JSON).
#
# CLI behavior:
#   • If you pass an input path: no GUI dialog is used.
#   • If you don't pass a path: we try to open a file dialog.
#   • If the dialog isn't available (headless) or is canceled: we print guidance and exit.
#   • --threads / --verbose flags override interactive prompts; if omitted, we still prompt.
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
from .settings import logger, ask_for_settings_cli, setup_logging

# Per-thread HTTP session (each worker keeps its own keep-alive pool)
_thread_local = threading.local()

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
    - If threads/verbose flags not provided: prompt interactively.
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
            "    demodifier your_input.csv\n"
            "Optionally add flags, e.g.:\n"
            "    demodifier your_input.csv --threads 8 --verbose"
        )
        return

    # 2) Resolve execution settings (prompt only for missing values)
    if num_threads is None and verbose is None:
        default_threads = max(1, (os.cpu_count() or 4))
        # Old behavior: ask both if neither provided
        num_threads, verbose = ask_for_settings_cli(default_threads, False)
    else:
        if num_threads is None:
            # Ask ONLY for processors
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

        if verbose is None:
            # Ask ONLY for verbose mode
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
                    verbose = False  # default if user just hits Enter
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
# Adds optional flags that *bypass* interactive prompts when provided.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process peptides from a CSV/TSV file.")
    parser.add_argument("input_csv", nargs="?", default=None, help="Path to input CSV/TSV file")
    parser.add_argument("--threads", type=int, help="Number of worker threads (overrides prompt)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging (overrides prompt)")
    args = parser.parse_args()

    # Pass flags through; if omitted, keep None so main() will prompt.
    main(
    input_csv=args.input_csv,
    num_threads=args.threads if args.threads is not None else None,
    verbose=args.verbose,
    )