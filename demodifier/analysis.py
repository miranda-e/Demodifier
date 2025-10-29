# analysis.py
# ---------------------------------------------------------
# The Demodifier generates peptide sequences and simulates potential sequence modifications (caused by deamidation, reamidation, pyroglu)
# to generate all possible Modification Induced Sequence Permutations (MISPs).
# The script then sends each MISP to the Unipept API to retrieve its lowest Common Ancestor (LCA).
# Finally, it outputs a results CSV, containing all MISPs and their LCAs,
# enabling the researcher to assess potentially incorrect peptide taxonomies.
# If you use this tool, please cite Evans (2025), 
# The Demodifier: a tool for screening modification-induced alternate peptide taxonomy in palaeoproteomics
# ---------------------------------------------------------

import re
import logging
from itertools import combinations

# Configure module logger (level/handlers managed by settings.setup_logging)
logger = logging.getLogger(__name__)

# Counts how many deamidation events are listed in the Modifications column
# Example: "2 Deamidated (NQ)" → returns 2
def extract_deamidation_count(modifications):
    if not modifications:
        return 0
    matches = re.findall(r'(?:(\d+)\s*)?(Deamidated|Deamidation)\s*\(NQ\)', modifications, re.IGNORECASE)
    return sum(int(num) if num else 1 for num, _ in matches)

# Count how many residues in the sequence can be deamidated (i.e. number of N or Q in peptide)
# Accounts for pyro-Glu shortening the sequence (i.e. ignore N-term Q or E if pyro-Glu mod was detected)
def count_deamidatable_residues(peptide, modifications):
    if not peptide:
        return 0
    pyro_glu = "pyro" in modifications.lower() if modifications else False
    sequence = peptide[1:] if pyro_glu and peptide[0] in {"Q", "E"} else peptide
    return sum(1 for aa in sequence if aa in {"N", "Q"})

# Decide if we need to care about exact deamidation position
# (only required if multiple MISPs give different LCAs)
def deamidation_position_required(modifications, peptide, identical_lcas):
    if not peptide:
        return "not required"
    num_deamid = extract_deamidation_count(modifications)
    if identical_lcas != "no" or num_deamid == 0:
        return "not required"
    nq_total = count_deamidatable_residues(peptide, modifications)
    return "required" if num_deamid < nq_total else "not required"

# Generate all possible deamidation-based permutations of the peptide
# in which up to the maximum number of deamidations detected are made (substitute N→D and Q→E)
# Ignores N-terminal Q or E if pyro-glu modification is detected (by inserting a placeholder amino acid, "X")
def generate_deamidation_permutations(peptide, max_substitutions, modifications, verbose=False):
    # Ensure modifications is a string, even if it's passed as a boolean
    modifications = str(modifications) if modifications else ""
    
    # Check for pyro-Glu modification
    pyro_glu = "pyro" in modifications.lower() if modifications else False
    
    # Pretend that the first residue (Q or E) is a placeholder amino acid 'X' for substitution 
    # to prevent inaccurate n-term subs (doesn't change actual peptide)
    peptide_for_permutation = peptide
    if pyro_glu and peptide[0] in {"Q", "E"}:
        peptide_for_permutation = "X" + peptide[1:]

    # Find all N and Q indices in the (modified X containing) peptide
    n_indices = [i for i, letter in enumerate(peptide_for_permutation) if letter == "N"]
    q_indices = [i for i, letter in enumerate(peptide_for_permutation) if letter == "Q"]
    
    if not n_indices and not q_indices:
        return [(peptide, 0, [])]  # No permutations if no N or Q in peptide
    
    permutations = []
    
    # Generate deamidation-based permutations (N -> D and Q -> E)
    for num_n in range(min(max_substitutions, len(n_indices)) + 1):
        for num_q in range(min(max_substitutions - num_n, len(q_indices)) + 1):
            for combo_n in combinations(n_indices, num_n):
                for combo_q in combinations(q_indices, num_q):
                    temp_peptide = list(peptide)  
                    
                    # Modify the positions with N -> D and Q -> E (on the full non-X-containing peptide)
                    modified_positions = list(combo_n) + list(combo_q)
                    for index in combo_n:
                        temp_peptide[index] = "D"  # Substitute N -> D
                    for index in combo_q:
                        temp_peptide[index] = "E"  # Substitute Q -> E
                    
                    # Store the permutation and the positions modified
                    permutations.append(("".join(temp_peptide), len(modified_positions), modified_positions))
    
    # Print all permutations if verbose flag is set
    if verbose:
        print(f"Deamidation permutations for peptide {peptide}:")
        for perm in permutations:
            print(perm[0])

    logger.debug(f"Generated {len(permutations)} deamidation-induced permutations.")
    return permutations

# Generate all reamidated-based permutations by substituting D→N and E→Q,
# except at positions previously substituted
# Ignores N-terminal Q or E if pyro-Glu modification is detected (by inserting a placeholder amino acid, "X")
def generate_reamidation_permutations(peptide, modified_positions, modifications=None, verbose=False):
    modifications = str(modifications) if modifications else ""

    # Check for pyro-Glu modification
    pyro_glu = "pyro" in modifications.lower() if modifications else False

    # If pyro-Glu modification is present, treat the first residue (Q or E) as a placeholder ('X') for permutation making
    # (doesn't alter actual peptide)
    peptide_for_permutation = peptide
    if pyro_glu and peptide[0] in {"Q", "E"}:
        peptide_for_permutation = "X" + peptide[1:]

    # Find all D and E indices in the (modified X containing) peptide, excluding positions already modified
    d_indices = [i for i, letter in enumerate(peptide_for_permutation) if letter == "D" and i not in modified_positions]
    e_indices = [i for i, letter in enumerate(peptide_for_permutation) if letter == "E" and i not in modified_positions]

    if not d_indices and not e_indices:
        return [peptide]

    permutations = []
    for num_d in range(len(d_indices) + 1):
        for num_e in range(len(e_indices) + 1):
            for combo_d in combinations(d_indices, num_d):
                for combo_e in combinations(e_indices, num_e):
                    temp_peptide = list(peptide)  # working with full (non "X"-containing) peptide
                    for index in combo_d:
                        temp_peptide[index] = "N"  # Substitute D -> N
                    for index in combo_e:
                        temp_peptide[index] = "Q"  # Substitute E -> Q
                    permutations.append("".join(temp_peptide))

    # Print permutations if verbose flag is set
    if verbose:
        print(f"Reamidation permutations for peptide {peptide}:")
        for perm in permutations:
            print(perm)

    return permutations

# Generate pyro-glu based permutations, if pyro-Glu conversions are specified in Modifications column
# I.e. if Q or E at the start of the sequence and Gln->pyro-Glu or Glu->pyro-Glu detected in search, 
# then simulate Q→E or E→Q respectively
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
