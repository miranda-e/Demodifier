# Updated demodifier with richer JSON output, logging, runtime tracking, and deamidation position checking flag

import argparse
import csv
import json
import re
import os
import time
import logging
import requests
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ---------------------------------------------
# Logging Configuration
# ---------------------------------------------
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)s] %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
# ---------------------------------------------
# Shared Session for HTTP Requests
# ---------------------------------------------
session = requests.Session()
retries = Retry(
    total=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    raise_on_status=False
)
session.mount("https://", HTTPAdapter(max_retries=retries))
# ---------------------------------------------
# API Call with Retry Logic
# ---------------------------------------------
def process_peptides(peptides, session=session):
    url = "https://api.unipept.ugent.be/api/v1/pept2lca"
    params = {"input[]": peptides, "equate_il": "true"}

    try:
        response = session.get(url, params=params, timeout=10)
        response.raise_for_status()
        results = response.json()
        lca_map = {result['peptide']: result.get('taxon_name', "no match") for result in results}
        return [lca_map.get(peptide, "no match") for peptide in peptides]
    except requests.exceptions.RequestException as e:
        logger.error(f"API request failed: {e}")
    except ValueError as e:
        logger.error(f"Failed to decode JSON from Unipept API: {e}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
    return ["no response"] * len(peptides)

# ---------------------------------------------
# Deamidation Utilities
# ---------------------------------------------
def extract_deamidation_count(modifications):
    if not modifications:
        return 0
    matches = re.findall(r'(?:(\d+)\s*)?(Deamidated|Deamidation)\s*\(NQ\)', modifications, re.IGNORECASE)
    count = 0
    for number, _ in matches:
        count += int(number) if number else 1
    return count

def count_deamidatable_residues(peptide, modifications):
    if not peptide:
        return 0
    pyro_glu = "pyro" in modifications.lower() if modifications else False
    sequence = peptide[1:] if pyro_glu and peptide[0] in {"Q", "E"} else peptide
    return sum(1 for aa in sequence if aa in {"N", "Q"})

def deamidation_position_required(modifications, peptide, identical_lcas):
    num_deamid = extract_deamidation_count(modifications)
    if identical_lcas != "no" or num_deamid == 0:
        return "not required"
    nq_total = count_deamidatable_residues(peptide, modifications)
    return "required" if num_deamid < nq_total else "not required"

# ---------------------------------------------
# Permutation Functions
# ---------------------------------------------
def generate_peptide_permutations(peptide, max_substitutions):
    n_indices = [i for i, letter in enumerate(peptide) if letter == "N"]
    q_indices = [i for i, letter in enumerate(peptide) if letter == "Q"]
    if not n_indices and not q_indices:
        return [(peptide, 0, [])]
    permutations = []
    for num_n in range(min(max_substitutions, len(n_indices)) + 1):
        for num_q in range(min(max_substitutions - num_n, len(q_indices)) + 1):
            for combo_n in combinations(n_indices, num_n):
                for combo_q in combinations(q_indices, num_q):
                    temp_peptide = list(peptide)
                    modified_positions = []
                    for index in combo_n:
                        temp_peptide[index] = "D"
                        modified_positions.append(index)
                    for index in combo_q:
                        temp_peptide[index] = "E"
                        modified_positions.append(index)
                    permutations.append(("".join(temp_peptide), len(modified_positions), modified_positions))
    logger.debug(f"Generated {len(permutations)} peptide permutations for {peptide}")
    return permutations

def generate_reamidation_permutations(peptide, modified_positions):
    d_indices = [i for i, letter in enumerate(peptide) if letter == "D" and i not in modified_positions]
    e_indices = [i for i, letter in enumerate(peptide) if letter == "E" and i not in modified_positions]
    if not d_indices and not e_indices:
        return [peptide]
    permutations = []
    for num_d in range(len(d_indices) + 1):
        for num_e in range(len(e_indices) + 1):
            for combo_d in combinations(d_indices, num_d):
                for combo_e in combinations(e_indices, num_e):
                    temp_peptide = list(peptide)
                    for index in combo_d:
                        temp_peptide[index] = "N"
                    for index in combo_e:
                        temp_peptide[index] = "Q"
                    permutations.append("".join(temp_peptide))
    logger.debug(f"Generated {len(permutations)} reamidation permutations for {peptide}")
    return permutations

def add_pyro_glu_permutations(peptide_permutations, modifications):
    pyro_permutations = set(peptide_permutations)
    original_peptide = peptide_permutations[0]
    modified_peptides = []
    modification_applied = False
    if modifications:
        if original_peptide.startswith("Q") and "Gln->pyro-Glu" in modifications:
            modified_peptides = [perm.replace("Q", "E", 1) for perm in peptide_permutations]
            modification_applied = True
        elif original_peptide.startswith("E") and "Glu->pyro-Glu" in modifications:
            modified_peptides = [perm.replace("E", "Q", 1) for perm in peptide_permutations]
            modification_applied = True
        for perm in modified_peptides:
            if perm != original_peptide:
                logger.debug(f"Added pyro-Glu permutation: {perm}")
        pyro_permutations.update(modified_peptides)
    if not modification_applied:
        logger.debug(f"No pyro-Glu modifications were applied to {original_peptide}")
    return list(pyro_permutations)

# ---------------------------------------------
# Row Processor
# ---------------------------------------------
def process_row(row):
    peptide = row.get('Sequence') or row.get('pep_seq')
    modifications = row.get('Modifications') or row.get('pep_var_mod')
    max_substitutions = extract_deamidation_count(modifications)
    peptide_options = generate_peptide_permutations(peptide, max_substitutions)
    final_permutations = []
    for perm, _, modified_positions in peptide_options:
        final_permutations.extend(generate_reamidation_permutations(perm, modified_positions))
    final_permutations = list(set(add_pyro_glu_permutations(final_permutations, modifications)))
    input_pep_lca = process_peptides([peptide])[0]
    return row, final_permutations, input_pep_lca

# ---------------------------------------------
# LCA Batch Processing
# ---------------------------------------------
def process_lca_in_batches(peptide_batches, max_workers=4):
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_peptides, batch) for batch in peptide_batches]
        for f in as_completed(futures):
            results.extend(f.result())
    return results

# ---------------------------------------------
# Main Program Execution
# ---------------------------------------------
def main(input_csv, num_threads, verbose=False):
    if verbose:
        logger.setLevel(logging.DEBUG)

    start_time = time.time()

    input_name = os.path.splitext(input_csv)[0]
    output_csv = f"{input_name}_results.csv"
    output_json = f"{input_name}_output.json"
    output_variant_lca_csv = f"{input_name}_permutations.csv"

    logger.info(f"Starting demodifier with {num_threads} threads...")

    with open(input_csv, 'r', encoding='utf-8-sig') as csvfile, \
         open(output_variant_lca_csv, 'w', newline='') as variantlcafile:

        reader = csv.DictReader(csvfile)
        variant_lca_writer = csv.writer(variantlcafile)
        variant_lca_writer.writerow(["Sequence", "Modifications", "Variant", "Variant_LCA"])

        fieldnames = [
            'Sequence', 'Modifications', 'input_pep_LCA', 'total_unique_permutations_(count)',
            'all_permutations_with_LCAs_(Count)', 'permutations_yeilding_LCAs',
            'all_permutation_LCAs', 'identical_LCAs', 'deamidation_position_checking'
        ]

        with open(output_csv, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()

            data = []
            peptide_to_lca = []
            rows_to_process = []

            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                futures = [executor.submit(process_row, row) for row in reader]
                for future in as_completed(futures):
                    row, final_permutations, input_pep_lca = future.result()
                    rows_to_process.append((row, final_permutations, input_pep_lca))
                    peptide_to_lca.extend(final_permutations)

            batch_size = 100
            peptide_batches = [peptide_to_lca[i:i + batch_size] for i in range(0, len(peptide_to_lca), batch_size)]
            lca_results = process_lca_in_batches(peptide_batches, max_workers=num_threads)

            lca_index = 0
            for row, final_permutations, input_pep_lca in rows_to_process:
                sequence_col = row.get('Sequence') or row.get('pep_seq')
                modifications_col = row.get('Modifications') or row.get('pep_var_mod')
                lca_options = []
                permutations_with_lcas = []
                permutation_json = []

                row_lcas = process_peptides(final_permutations)
                for variant, lca in zip(final_permutations, row_lcas):
                    variant_lca_writer.writerow([sequence_col, modifications_col, variant, lca])
                    permutation_json.append({"variant": variant, "lca": lca})
                    if lca != "no match":
                        lca_options.append(lca)
                    permutations_with_lcas.append(variant)

                identical_lcas = "NA" if len(lca_options) <= 1 else ("yes" if len(set(lca_options)) == 1 else "no")
                deamid_check = deamidation_position_required(modifications_col, sequence_col, identical_lcas)

                output_row = {
                    'Sequence': sequence_col,
                    'Modifications': modifications_col,
                    'input_pep_LCA': input_pep_lca,
                    'total_unique_permutations_(count)': len(final_permutations),
                    'all_permutations_with_LCAs_(Count)': len(permutations_with_lcas),
                    'permutations_yeilding_LCAs': ";".join(permutations_with_lcas),
                    'all_permutation_LCAs': ";".join(lca_options),
                    'identical_LCAs': identical_lcas,
                    'deamidation_position_checking': deamid_check
                }

                writer.writerow(output_row)

                data.append({
                    "sequence": sequence_col,
                    "modifications": modifications_col,
                    "input_pep_LCA": input_pep_lca,
                    "permutations": permutation_json,
                    "summary": {
                        "total_permutations": len(final_permutations),
                        "matched_permutations": len(permutations_with_lcas),
                        "unique_LCAs": list(set(lca_options)),
                        "identical_LCAs": identical_lcas == "yes",
                        "deamidation_position_checking": deamid_check
                    }
                })

    with open(output_json, 'w', encoding='utf-8') as jsonfile:
        json.dump(data, jsonfile, ensure_ascii=False, indent=4)

    elapsed = time.time() - start_time
    logger.info(f"Script completed in {elapsed:.2f} seconds using {num_threads} threads.")

# ---------------------------------------------
# CLI
# ---------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process peptides from a CSV file.")
    parser.add_argument("input_csv", help="Input CSV file containing peptide sequences")
    parser.add_argument("--num-processors", type=int, default=4, help="Number of threads to use for processing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug output")
    args = parser.parse_args()
    main(args.input_csv, args.num_processors, args.verbose)
