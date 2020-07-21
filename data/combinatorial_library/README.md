# KinFragLib: Combinatorial library 

This folder is meant for the metadata and properties of the KinFragLib combinatorial library, which is based on the KinFragLib fragment library at https://github.com/volkamerlab/KinFragLib. This dataset is used for the analysis of the combinatorial library.

Since this dataset contains large files, we provide it outside this repository at https://osf.io/cy7tr/.
In order to run the analysis notebooks, please download this dataset to this folder (as also instructed in the respective notebooks). 

## Raw data (download from https://osf.io/cy7tr/)

- `combinatorial_library.json`: Full combinatorial library, please refer to `notebooks/4_1_combinatorial_library_data_preparation.ipynb` at https://github.com/volkamerlab/KinFragLib for detailed information about this data format
- `combinatorial_library_deduplicated.json`: Deduplicated combinatorial library (based on InChIs)
- `chembl_standardized_inchi.csv`: Standardized ChEMBL 25 molecules in the form of InChI strings.

## Processed data (download from https://osf.io/cy7tr/)

Data extracted from `combinatorial_library_deduplicated.json`, performed in `notebooks/4_1_combinatorial_library_data_preparation.ipynb` at https://github.com/volkamerlab/KinFragLib.

- `metadata.json`: Bond and fragment IDs for each recombined ligand
- `n_atoms.csv`: Number of atoms for each recombined ligand
- `ro5.csv`: Number of ligands that fulfill Lipinski's rule of five (Ro5) and its individual criteria; number of ligands in total
- `subpockets.csv`: Number of ligands per subpocket combination
- `original_exact.json`: Ligands with exact matches in original ligands, i.e. KLIFS ligands that were used for the fragmentation
- `original_substructure.json`: Ligands with substructure matches in original ligands, i.e. KLIFS ligands that were used for the fragmentation
- `chembl_exact.json`: Ligands with exact matches in ChEMBL
- `chembl_most_similar.json`: Most similar ligand in ChEMBL for each recombined ligand 


