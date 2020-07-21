# Notebooks

Overview on notebook content:

## 1. Quick start

### `1_1_quick_start.ipynb`

This notebook contains an introduction on how to load and use the KinFragLib fragment library.

## 2. Full fragment library

### `2_1_fragment_analysis_original_ligands.ipynb`

Based on the full fragment library, get all ligands from which these fragments originated from ("original ligands" from KLIFS dataset).

### `2_2_fragment_analysis_statistics.ipynb`

This notebook contains the code that was used to calculate most of the statistics as well as to generate the respective plots shown in the manuscript, for instance:

- Ligand occupancy across subpockets
- Subpocket connectivity across subpockets
- Fragment occurrence per subpocket
- Fragment properties per subpocket
- Fragment similarity per subpocket
- Fragment promiscuity 

### `2_3_fragment_analysis_most_common_fragments.ipynb`

This notebook contains the code that was used to analyze the most common fragments per subpocket as well as to generate the respective plots shown in the manuscript.

- Find 50 most common fragments in each subpocket (if multiple fragments with same count are at cutoff, include all fragments).
- Cluster these fragments using Butina clustering.
- Draw 50 most common fragments per subpocket sorted by descending cluster size.

### `2_4_fragment_analysis_literature.ipynb`

[WIP] This notebook aims to load kinase inhibitor fragment data from literature and check if these fragments are included in our fragment library.

## 3. Reduced fragment library

### `3_1_fragment_library_reduced.ipynb`

The fragment library resulting from the KinFragLib fragmentation procedure comprises of about 3000 fragments. Ultimately, we want to demonstrate how this library can be used for recombining ligands. Before this can be done, we need to address two considerations:

1. Remove all fragments that are not useful in a recombination, i.e. duplicates, fragments in pool X, fragments without dummy atoms, and fragments with dummy atoms only connecting to pool X. Also remove all AP fragments that show no hydrogen bond donors and acceptors (not hinge-like).
2. Select a diverse set of fragments (per subpocket) for recombination to (i) save computational cost and (ii) avoid recombination of highly similar fragments.

## 4. Combinatorial library

### `4_1_combinatorial_library_data.ipynb`

The aim of this notebook is to extract information from the combinatorial library (`json` file) about e.g. ligand sizes, Lipinski's rule of five compliance, and matches in ChEMBL and KLIFS. Since the `json` file holds mulitple millions of ligands, we do this data processing once here at the beginning and save the results to separate files which will be used for analysis/visualization in the following notebooks.

### `4_2_combinatorial_library_properties.ipynb`

In this notebook we want to analyze properties of the combinatorial library, such as the ligand size and Lipinski's rule of five criteria.

### `4_3_combinatorial_library_comparison_klifs.ipynb`

In this notebook we want to compare the combinatorial library to the original KLIFS ligands, i.e. the ligands from which the fragment library originates from. We consider exact and substructure matches.

### `4_4_combinatorial_library_comparison_chembl.ipynb`

In this notebook we want to compare the combinatorial library to the ChEMBL 25 dataset in order to find exact matches and the most similar ChEMBL ligand per recombined ligand.
