# Filtered fragment library

The (full) fragment library resulting from the KinFragLib fragmentation procedure comprises of 9505 fragments, which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases (see `data/fragment_library/`).

In order to prepare a library with fragments tailored for recombination, we offer heare a filtered fragment library (2383 fragments) based on the following filters:

1. Remove pool X
2. Deduplicate fragment library (per subpocket)
3. Remove fragments without dummy atoms (unfragmented ligands)
4. Remove all fragments that connect only to pool X
5. Keep "Rule of Three (Ro3)" compliant fragments (fragment-likeness)
6. Filter AP subpocket fragments (typical hinge-like)

Please refer to the notebook `notebooks/3_1_fragment_library_reduced.ipynb` to check how the data was generated and/or use it as a starting point for customized protocols.