import argparse
import csv
import json
import requests
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import os

# Batch processing API function
def process_peptides(peptides):
    url = "https://api.unipept.ugent.be/api/v1/pept2lca"
    params = {"input[]": peptides, "equate_il": "true"}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        results = response.json()
        lca_map = {result['peptide']: result.get('taxon_name', "no match") for result in results}
        return [lca_map.get(peptide, "no match") for peptide in peptides]
    else:
        return ["no match"] * len(peptides)

# Generate peptide permutations (Phase 1 - N to D and Q to E)
def generate_peptide_permutations(peptide, max_substitutions):
    n_indices = [i for i, letter in enumerate(peptide) if letter == "N"]
    q_indices = [i for i, letter in enumerate(peptide) if letter == "Q"]

    if not n_indices and not q_indices:
        return [(peptide, 0)]  # Return original peptide with 0 substitutions

    permutations = []
    
    for num_n in range(min(max_substitutions, len(n_indices)) + 1):
        for num_q in range(min(max_substitutions - num_n, len(q_indices)) + 1):
            for combo_n in combinations(n_indices, num_n):
                for combo_q in combinations(q_indices, num_q):
                    temp_peptide = list(peptide)
                    substitution_count = 0

                    for index in combo_n:
                        temp_peptide[index] = "D"
                        substitution_count += 1

                    for index in combo_q:
                        temp_peptide[index] = "E"
                        substitution_count += 1

                    permutations.append(("".join(temp_peptide), substitution_count))

    return permutations

# Generate reamidation permutations (Phase 2 - D to N and E to Q)
def generate_reamidation_permutations(peptide):
    d_indices = [i for i, letter in enumerate(peptide) if letter == "D"]
    e_indices = [i for i, letter in enumerate(peptide) if letter == "E"]

    if not d_indices and not e_indices:
        return [peptide]

    permutations = []

    for num_d in range(len(d_indices) + 1):  # Unlimited substitutions for reamidation
        for num_e in range(len(e_indices) + 1):  # Unlimited substitutions for reamidation
            for combo_d in combinations(d_indices, num_d):
                for combo_e in combinations(e_indices, num_e):
                    temp_peptide = list(peptide)

                    for index in combo_d:
                        temp_peptide[index] = "N"

                    for index in combo_e:
                        temp_peptide[index] = "Q"

                    permutations.append("".join(temp_peptide))

    return permutations

# Additional permutation step for pyro-Glu conversion
def add_pyro_glu_permutations(peptide_permutations, modifications):
    pyro_permutations = set(peptide_permutations)  # Use set to avoid duplicates

    if modifications:
        if peptide_permutations[0].startswith("Q") and "Gln->pyro-Glu" in modifications:
            pyro_permutations.update([perm.replace("Q", "E", 1) for perm in peptide_permutations])
        elif peptide_permutations[0].startswith("E") and "Glu->pyro-Glu" in modifications:
            pyro_permutations.update([perm.replace("E", "Q", 1) for perm in peptide_permutations])

    return list(pyro_permutations)

# Function to extract the deamidation count from the modifications string
def extract_deamidation_count(modifications):
    if not modifications:
        return 0
    deamidation_match = re.search(r'(\d+)?\s*(?:Deamidated|Deamidation)\s*\(NQ\)', modifications)
    if deamidation_match:
        return int(deamidation_match.group(1)) if deamidation_match.group(1) else 1
    return 0

# Function to process a single row
def process_row(row):
    peptide = row.get('Sequence') or row.get('pep_seq')
    modifications = row.get('Modifications') or row.get('pep_var_mod')

    max_substitutions = extract_deamidation_count(modifications)
    peptide_options = generate_peptide_permutations(peptide, max_substitutions)
    final_permutations = []

    for perm, count in peptide_options:
        reamidation_permutations = generate_reamidation_permutations(perm)
        final_permutations.extend(reamidation_permutations)

    final_permutations = list(set(add_pyro_glu_permutations(final_permutations, modifications)))

    # Get input peptide LCA
    input_pep_lca = process_peptides([peptide])[0]  # Retrieve LCA for the original peptide

    return row, final_permutations, input_pep_lca

# Function to handle LCA processing in batches
def process_lca_in_batches(peptide_batches):
    batch_results = []
    for batch in peptide_batches:
        batch_results.extend(process_peptides(batch))
    return batch_results

# Main function to run the program
def main(input_csv):
    input_name = os.path.splitext(input_csv)[0]
    output_csv = f"{input_name}_results.csv"
    output_json = f"{input_name}_output.json"
    output_variant_lca_csv = f"{input_name}_permutations.csv"
    
    with open(input_csv, 'r', encoding='utf-8-sig') as csvfile, \
         open(output_variant_lca_csv, 'w', newline='') as variantlcafile:
        
        reader = csv.DictReader(csvfile)
        variant_lca_writer = csv.writer(variantlcafile)
        variant_lca_writer.writerow(["Sequence", "Modifications", "Variant", "Variant_LCA"])
        
        fieldnames = [
            'Sequence', 'Modifications', 'input_pep_LCA', 'total_unique_permutations_(count)', 
            'all_permutations_with_LCAs_(Count)', 'Permutations Yielding LCAs', 
            'all_permutation_LCAs', 'Identical_LCAs' 
        ]

        with open(output_csv, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            
            data = []
            peptide_to_lca = []
            rows_to_process = []

            # Step 1: Parallel row processing
            with ThreadPoolExecutor() as executor:
                futures = []
                for row in reader:
                    futures.append(executor.submit(process_row, row))

                for future in as_completed(futures):
                    row, final_permutations, input_pep_lca = future.result()
                    rows_to_process.append((row, final_permutations, input_pep_lca))

                    peptide_to_lca.extend(final_permutations)

            # Step 2: Batch process peptides for LCA
            batch_size = 100  
            peptide_batches = [peptide_to_lca[i:i + batch_size] for i in range(0, len(peptide_to_lca), batch_size)]
            lca_results = process_lca_in_batches(peptide_batches)

            # Step 3: Write final outputs
            lca_index = 0
            for row, final_permutations, input_pep_lca in rows_to_process:
                lca_options = []
                permutations_with_lcas = []
                
                for variant in final_permutations:
                    lca = lca_results[lca_index]
                    lca_index += 1

                    sequence_col = row.get('Sequence') or row.get('pep_seq')
                    modifications_col = row.get('Modifications') or row.get('pep_var_mod')

                    variant_lca_writer.writerow([sequence_col, modifications_col, variant, lca])

                    if lca != "no match":
                        lca_options.append(lca)
                        permutations_with_lcas.append(variant)

                row['total_unique_permutations_(count)'] = len(final_permutations)
                row['all_permutations_with_LCAs_(Count)'] = len(permutations_with_lcas)
                row['Permutations Yielding LCAs'] = ";".join(permutations_with_lcas)
                row['all_permutation_LCAs'] = ";".join(lca_options)
                row['input_pep_LCA'] = input_pep_lca

                if len(lca_options) <= 1:
                    row['Identical_LCAs'] = "NA"
                elif len(set(lca_options)) == 1:
                    row['Identical_LCAs'] = "yes"
                else:
                    row['Identical_LCAs'] = "no"

                output_row = {
                    'Sequence': sequence_col,
                    'Modifications': modifications_col,
                    'input_pep_LCA': row['input_pep_LCA'],  # Include the new column
                    'total_unique_permutations_(count)': row['total_unique_permutations_(count)'],
                    'all_permutations_with_LCAs_(Count)': row['all_permutations_with_LCAs_(Count)'],
                    'Permutations Yielding LCAs': row['Permutations Yielding LCAs'],
                    'all_permutation_LCAs': row['all_permutation_LCAs'],
                    'Identical_LCAs': row['Identical_LCAs']  
                }
                writer.writerow(output_row)

    with open(output_json, 'w', encoding='utf-8') as jsonfile:
        json.dump(data, jsonfile, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process peptides from a CSV file.")
    parser.add_argument("input_csv", help="Input CSV file containing peptide sequences")
    args = parser.parse_args()
    main(args.input_csv)
