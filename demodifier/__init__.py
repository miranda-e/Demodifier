# __init__.py
# ---------------------------------------------------------
# The Demodifier simulates potential sequence modifications (caused by deamidation and pyroglu formation)
# to generate all possible Modification Induced Sequence Permutations (MISPs).
# The script then sends each MISP to the Unipept API to retrieve its lowest Common Ancestor (LCA).
# Finally, it outputs a results CSV, containing all MISPs and their LCAs,
# enabling the researcher to assess potentially incorrect peptide taxonomies.
# If you use this tool, please cite Evans (2025), 
# The Demodifier: a tool for screening modification-induced alternate peptide taxonomy in palaeoproteomics
# ---------------------------------------------------------

# Intentionally empty to avoid importing submodules at package import time.
