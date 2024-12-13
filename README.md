# Demodifier
The Demodifier generates potential modification-induced alternate peptide sequences and their lowest common taxonomic ancestor given a csv containing a list of peptide sequences and their post translational modifications, to avoid spurious taxonomic detections in ancient protein studies.

The Demodifier is a python script which screens for possible modification-induced sequence permutations (currently supporting deamidation of N and Q and pyroglutamic acid formation at N-terminus E and Q), and calls the Unipept pept2lca API (Mesuere et al. 2016) to assign each peptide permutation its lowest common ancestor, allowing the researcher to scrutinise possibly inaccurate peptide taxonomies in ancient protein studies.

Given a csv containing a list of peptide sequences in one column, and their accompanying modifications in Mascot or Maxquant format in another, The Demodifier will simulate all possible modification induced sequence permutations, and call the Unipept pept2lca api to assign taxonomy to each. This enables the researcher to scrutinise all possible LCAs for a given peptide, avoiding spurious taxonomic detections caused by modifications.

To run The Demodifier, save the script and your input csv in the same directory, and then run the script using command:

Python3 De-modifier.py your_input_file.csv

Please refer to Evans et al. (2024) for further information.

If you use this tool, please cite:

Evans (2024): The Demodifier: a tool for screening modification-induced alternate peptide taxonomy in palaeoproteomics (preprint)

Bart Mesuere, Toon Willems, Felix Van der Jeugt, Bart Devreese, Peter Vandamme, Peter Dawyndt, Unipept web services for metaproteomics analysis, Bioinformatics, Volume 32, Issue 11, June 2016, Pages 1746–1748, https://doi.org/10.1093/bioinformatics/btw039
