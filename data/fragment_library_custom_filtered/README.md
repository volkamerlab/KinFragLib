# Custom filtered fragment library

The (full) fragment library resulting from the KinFragLib fragmentation procedure comprises of 7486 fragments, which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases (see `data/fragment_library/`).

To reduce the fragment library size and enable the recombination avoiding the recombinatorial explosion we offer the possibility to filter the fragment library by customizable filtering steps, namely:

1. Pre-filtering (Remove pool X, deduplicate, remove unfragmented fragments, remove fragments only connecting to pool X and fragments in Pool X) \[not optional\]
2. Filter for unwanted substructures (PAINS and Brenk et al.) \[optional\]
3. Filter for drug likeness (Ro3 and QED) \[optional\]
4. Filter for synthesizability (Buyable building blocks and SYBA) \[optional\]
5. filter for pairwise retrosynthesizability (using ASKCOS) \[optional\]

Please refer to the notebook `notebooks/custom_kinfraglib/2_1_custom_filters_pipeline.ipynb` to check how the data was generated and/or to generate your own custom fragment library (de-)activating filters and modifying the filters parameters.


- `custom_filter_results.csv`: file containing the filtering results from the analysis applied, including the fragments SMILES and the fragments subpocket as indices and the values and boolean columns for each fragment (from the pre-filtered library) retrieved by the filtering steps.