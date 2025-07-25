# The Demodifier, Miranda Evans 2025
# ---------------------------------------------------------
# The Demodifier generates peptide sequences and simulatulates potential sequence modifications (caused by deamidation, reamidation, pyroglu)
# to generate all possible Modification Induced Sequence Permutations (MISPs).
# The script then sends each MISP to the Unipept API to retrieve its lowest Common Ancestor (LCA).
# Finally, it outputs a results CSV, containing all MISPs and their LCAs,
# enabling the researcher to assess potentially incorrect peptide taxonomies.
# If you use this tool, please cite Evans (2025), 
# The Demodifier: a tool for screening modification-induced alternate peptide taxonomy in palaeoproteomics
# ---------------------------------------------------------

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

# Configure logging to show status info
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)s] %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

# Setup a resilient web session (used for API calls)
# Handles retries on failure and avoids re-opening connections
session = requests.Session()
retries = Retry(
    total=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    raise_on_status=False
)
session.mount("https://", HTTPAdapter(max_retries=retries))

# Query the Unipept API with a list of peptide strings
# Returns their LCA
def process_peptides(peptides, session=session):
    url = "https://api.unipept.ugent.be/api/v1/pept2lca"
    params = {"input[]": peptides, "equate_il": "true"}
    
    # Debug: Log the peptides being queried
    logger.debug(f"Querying Unipept API with peptides: {peptides}")
    
    try:
        response = session.get(url, params=params, timeout=10)
        response.raise_for_status()
        results = response.json()
        
        # Debug: Log the raw results from the API
        logger.debug(f"API response: {results}")
        
        lca_map = {result['peptide']: result.get('taxon_name', "no match") for result in results}
        return [lca_map.get(peptide, "no match") for peptide in peptides]
    
    except Exception as e:
        logger.error(f"API request failed: {e}")
        return ["no response"] * len(peptides)

# Counts how many deamidation events are listed in the Modifications column
# Example: "2 Deamidated (NQ)" → returns 2
def extract_deamidation_count(modifications):
    if not modifications:
        return 0
    matches = re.findall(r'(?:(\d+)\s*)?(Deamidated|Deamidation)\s*\(NQ\)', modifications, re.IGNORECASE)
    return sum(int(num) if num else 1 for num, _ in matches)

# Count how many residues in the sequence can be deamidated (i.e. number of N or Q in peptide)
# Accounts for pyro-Glu shortening the sequence (i.e. ignore N-term Q or E if pyroGlu mod detected)
def count_deamidatable_residues(peptide, modifications):
    if not peptide:
        return 0
    pyro_glu = "pyro" in modifications.lower() if modifications else False
    sequence = peptide[1:] if pyro_glu and peptide[0] in {"Q", "E"} else peptide
    return sum(1 for aa in sequence if aa in {"N", "Q"})

# Decide if we need to care about exact deamidation position
# (only required if multiple variants give different LCAs)
def deamidation_position_required(modifications, peptide, identical_lcas):
    if not peptide:
        return "not required"
    num_deamid = extract_deamidation_count(modifications)
    if identical_lcas != "no" or num_deamid == 0:
        return "not required"
    nq_total = count_deamidatable_residues(peptide, modifications)
    return "required" if num_deamid < nq_total else "not required"

# Generate all possible deamidation-based permutations of the peptide
# where up to `max_substitutions` N→D and Q→E changes are made
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
                    modified_positions = list(combo_n) + list(combo_q)
                    for index in combo_n:
                        temp_peptide[index] = "D"
                    for index in combo_q:
                        temp_peptide[index] = "E"
                    permutations.append(("".join(temp_peptide), len(modified_positions), modified_positions))
    
    # Debug: Log number of permutations generated
    logger.debug(f"Generated {len(permutations)} deamidation-induced permutations.")
    return permutations

# Generate all reamidated-based permutations by reverting D→N and E→Q,
# except at positions previously substituted
def generate_reamidation_permutations(peptide, modified_positions, modifications=None):
    pyro_glu = "pyro" in modifications.lower() if modifications else False

    d_indices = [
        i for i, letter in enumerate(peptide)
        if letter == "D" and i not in modified_positions
    ]

    e_indices = [
        i for i, letter in enumerate(peptide)
        if letter == "E"
        and i not in modified_positions
        and not (i == 0 and pyro_glu)  # Skip N-terminal E if pyro-Glu is detected
    ]

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


    return permutations

# Generate pyro-glu based permutations, if pyro-Glu conversions are specified in Modifications column
# I.e. Q→E or E→Q at the start of the sequence and Gln->pyro-Glu or Glu->pyro-Glu detected in search
def add_pyro_glu_permutations(peptide_permutations, modifications):
    pyro_permutations = set(peptide_permutations)
    original_peptide = peptide_permutations[0]
    modified_peptides = []
    if modifications:
        if original_peptide.startswith("Q") and "Gln->pyro-Glu" in modifications:
            modified_peptides = [perm.replace("Q", "E", 1) for perm in peptide_permutations]
        elif original_peptide.startswith("E") and "Glu->pyro-Glu" in modifications:
            modified_peptides = [perm.replace("E", "Q", 1) for perm in peptide_permutations]
        pyro_permutations.update(modified_peptides)
    return list(pyro_permutations)

# Look up LCAs for a list of peptide variants in chunks of 100
def get_lcas_for_permutations(permutations, session):
    chunk_size = 100
    lca_map = {}
    for i in range(0, len(permutations), chunk_size):
        chunk = permutations[i:i + chunk_size]
        chunk_lcas = process_peptides(chunk, session=session)
        for pep, lca in zip(chunk, chunk_lcas):
            lca_map[pep] = lca
    return [lca_map.get(p, "no match") for p in permutations]

# Process one row (peptide + modifications) and return all variants + LCA data
def process_row(row, session):
    peptide = (row.get('Sequence') or row.get('pep_seq') or '').strip()
    modifications = (row.get('Modifications') or row.get('pep_var_mod') or '').strip()
    
    # Debug: Log peptide and modifications being processed
    logger.debug(f"Processing row: Peptide = {peptide}, Modifications = {modifications}")
    
    if not peptide:
        return None
    
    max_substitutions = extract_deamidation_count(modifications)
    
    # Debug: Log deamidation count
    logger.debug(f"Deamidation count for {peptide}: {max_substitutions}")
    
    peptide_options = generate_peptide_permutations(peptide, max_substitutions)
    
    # Debug: Log number of permutations generated
    logger.debug(f"Generated total {len(peptide_options)} MISPs for {peptide}")
    
    final_permutations = []
    for perm, _, modified_positions in peptide_options:
        final_permutations.extend(generate_reamidation_permutations(perm, modified_positions, modifications))
    
    # Debug: Log number of reamidation permutations
    logger.debug(f"Generated {len(final_permutations)} reamidation-induced permutations for {peptide}")
    
    final_permutations = list(set(add_pyro_glu_permutations(final_permutations, modifications)))
    
    # Debug: Log pyro-Glu permutations
    logger.debug(f"Generated {len(final_permutations)} pyro-Glu-induced permutations for {peptide}")
    
    input_pep_lca = process_peptides([peptide], session=session)[0]
    
    # Debug: Log LCA for the input peptide
    logger.debug(f"Input peptide LCA for {peptide}: {input_pep_lca}")
    
    row_lcas = get_lcas_for_permutations(final_permutations, session)
    
    return row, final_permutations, row_lcas, input_pep_lca

# Main script execution
# Loads input CSV, processes rows, writes outputs
def main(input_csv, num_threads, verbose=False):
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    start_time = time.time()
    input_name = os.path.splitext(input_csv)[0]
    output_csv = f"{input_name}_results.csv"
    output_json = f"{input_name}_output.json"
    output_variant_lca_csv = f"{input_name}_permutations.csv"
    
    # Debug: Log the script start
    logger.debug(f"Starting demodifier with {num_threads} threads on {input_csv}...")
    
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

            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                futures = [executor.submit(process_row, row, session) for row in reader]
                for future in as_completed(futures):
                    result = future.result()
                    if result is None:
                        continue
                    
                    row, final_permutations, row_lcas, input_pep_lca = result
                    
                    # Debug: Log progress on processing each row
                    logger.debug(f"Processing row: {row['Sequence']}")

                    sequence_col = (row.get('Sequence') or row.get('pep_seq') or '').strip()
                    modifications_col = (row.get('Modifications') or row.get('pep_var_mod') or '').strip()
                    lca_options = []
                    permutations_with_lcas = []
                    permutation_json = []

                    for variant, lca in zip(final_permutations, row_lcas):
                        variant_lca_writer.writerow([sequence_col, modifications_col, variant, lca])
                        permutation_json.append({"variant": variant, "lca": lca})
                        if lca != "no match":
                            lca_options.append(lca)
                            permutations_with_lcas.append(variant)

                    unique_lcas = set(lca_options)
                    identical_lcas = "NA" if len(lca_options) <= 1 else ("yes" if len(unique_lcas) == 1 else "no")
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
                            "unique_LCAs": list(unique_lcas),
                            "identical_LCAs": identical_lcas == "yes",
                            "deamidation_position_checking": deamid_check
                        }
                    })

    with open(output_json, 'w', encoding='utf-8') as jsonfile:
        json.dump(data, jsonfile, ensure_ascii=False, indent=4)

    elapsed = time.time() - start_time
    logger.debug(f"Script completed in {elapsed:.2f} seconds using {num_threads} threads.")

# Command-line interface: handles argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process peptides from a CSV file.")
    parser.add_argument("input_csv", help="Input CSV file containing peptide sequences")
    parser.add_argument("--num-processors", type=int, default=4, help="Number of threads to use for processing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug output")
    args = parser.parse_args()
    main(args.input_csv, args.num_processors, args.verbose)
