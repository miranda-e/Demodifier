# writers.py
# ---------------------------------------------------------
# Writes the Demodifier output files:
#   • Summary CSV: one row per input peptide
#   • Permutations CSV: one row per permutation (variant)
#   • JSON results: all information from summary and permutations CSVs, but in nested format
# No analysis or API logic lives here.
# ---------------------------------------------------------

import csv
import json
from .settings import logger
from .analysis import deamidation_position_required


def write_results(results, output_csv: str, output_json: str, output_variant_lca_csv: str):
    """
    Saves all results to CSV and JSON files.

    Parameters
    ----------
    results : list of tuples
        Each tuple contains:
          (1) row_dict             :original peptide data
          (2) final_permutations   :all generated variants
          (3) row_lcas             :LCA (lowest common ancestor) for each permutation (variant)
          (4) input_pep_lca        :LCA for the original input peptide
    output_csv : str
        File path for the summary CSV
    output_json : str
        File path for the JSON results
    output_variant_lca_csv : str
        File path for the permutations (variant) CSV
    """

    # Define the column headers for the summary CSV file
    fieldnames = [
        'Sequence',
        'Modifications',
        'input_pep_LCA',
        'total_unique_permutations_(count)',
        'all_permutations_with_LCAs_(count)',
        'permutations_yielding_LCAs',
        'all_permutation_LCAs',
        'identical_LCAs',
        'deamidation_position_checking'
    ]

    data = []   # This will hold all data for writing to the JSON file later
    # Open two output files at once:
    # - One for detailed variant results
    # - One for the overall summary
    with open(output_variant_lca_csv, 'w', newline='') as variantlcafile, \
         open(output_csv, 'w', newline='') as outfile:

        # Write the detailed per-variant CSV
        variant_lca_writer = csv.writer(variantlcafile)
        variant_lca_writer.writerow(["Sequence", "Modifications", "Variant", "Variant_LCA"])

        # Write the summary CSV
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()  # Write the column titles first
        # Loop through each peptide’s data in the results list
        for row, final_permutations, row_lcas, input_pep_lca in results:
            # Log progress for debugging or tracking
            logger.debug(f"Writing outputs for: {row.get('Sequence') or row.get('pep_seq')}")
            # Get peptide sequence and modification info, handling alternative key names
            sequence_col = (row.get('Sequence') or row.get('pep_seq') or '').strip()
            modifications_col = (row.get('Modifications') or row.get('pep_var_mod') or '').strip()
            # Prepare empty lists to collect data for this peptide
            lca_options = []
            permutations_with_lcas = []
            permutation_json = []
            # Go through each permutation (variant) and its corresponding LCA
            for variant, lca in zip(final_permutations, row_lcas):
                # Write a line in the variant-level CSV
                variant_lca_writer.writerow([sequence_col, modifications_col, variant, lca])

                # Add the pair to the JSON structure
                permutation_json.append({"variant": variant, "lca": lca})

                # Only count LCAs that are valid (not marked "no match")
                if lca != "no match":
                    lca_options.append(lca)
                    permutations_with_lcas.append(variant)
            # Figure out whether all LCAs are the same or different
            unique_lcas = set(lca_options)
            if len(lca_options) <= 1:
                identical_lcas = "NA"    # Not enough to compare (none or only one)
            elif len(unique_lcas) == 1:
                identical_lcas = "yes"   # All LCAs are identical
            else:
                identical_lcas = "no"    # There are different LCAs

            # Run the deamidation-position check function
            deamid_check = deamidation_position_required(modifications_col, sequence_col, identical_lcas)

            # Write one line to the summary CSV file
            writer.writerow({
                'Sequence': sequence_col,
                'Modifications': modifications_col,
                'input_pep_LCA': input_pep_lca,
                'total_unique_permutations_(count)': len(final_permutations),
                'all_permutations_with_LCAs_(count)': len(permutations_with_lcas),
                'permutations_yielding_LCAs': ";".join(permutations_with_lcas),
                'all_permutation_LCAs': ";".join(lca_options),
                'identical_LCAs': identical_lcas,
                'deamidation_position_checking': deamid_check
            })
            # Also prepare the same information for the JSON file
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

    # Once all peptides are processed, dump everything into the JSON file.
    with open(output_json, 'w', encoding='utf-8') as jsonfile:
        json.dump(data, jsonfile, ensure_ascii=False, indent=4)

    # Log file locations for verification/debugging
    logger.debug(
        f"Results written:\n"
        f"  • Summary CSV: {output_csv}\n"
        f"  • Variant LCA CSV: {output_variant_lca_csv}\n"
        f"  • JSON results: {output_json}"
    )