# writers.py
# ---------------------------------------------------------
# Writes the Demodifier outputs:
#   • Summary CSV (one row per input peptide)
#   • Variant LCA CSV (one row per permutation)
#   • JSON summary (machine-readable)
# No analysis or API logic lives here—just serialization.
# ---------------------------------------------------------

import csv
import json
from .settings import logger
from .analysis import deamidation_position_required


def write_results(results, output_csv: str, output_json: str, output_variant_lca_csv: str):
    """
    Serialize in-memory results to disk.

    Parameters
    ----------
    results : list[tuple]
        Each item: (row_dict, final_permutations, row_lcas, input_pep_lca)
    output_csv : str
        Path to summary CSV.
    output_json : str
        Path to JSON summary.
    output_variant_lca_csv : str
        Path to per-variant CSV.
    """
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

    data = []

    with open(output_variant_lca_csv, 'w', newline='') as variantlcafile, \
         open(output_csv, 'w', newline='') as outfile:

        # Detailed per-variant CSV
        variant_lca_writer = csv.writer(variantlcafile)
        variant_lca_writer.writerow(["Sequence", "Modifications", "Variant", "Variant_LCA"])

        # Summary CSV
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row, final_permutations, row_lcas, input_pep_lca in results:
            logger.debug(f"Writing outputs for: {row.get('Sequence') or row.get('pep_seq')}")

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
            if len(lca_options) <= 1:
                identical_lcas = "NA"
            elif len(unique_lcas) == 1:
                identical_lcas = "yes"
            else:
                identical_lcas = "no"

            deamid_check = deamidation_position_required(modifications_col, sequence_col, identical_lcas)

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

    logger.debug(
        f"Results written:\n"
        f"  • Summary CSV: {output_csv}\n"
        f"  • Variant LCA CSV: {output_variant_lca_csv}\n"
        f"  • JSON summary: {output_json}"
    )
