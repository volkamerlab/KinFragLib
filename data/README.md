# Data

**Note**: Our fragmentation library is currently based on KLIFS downloaded on *06.12.2023*. 

Overview of data content:

- `fragment_library/`: Full fragment library resulting from the KinFragLib fragmentation procedure comprises about 3000 fragments, which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases.
- `fragment_library_filtered/`: Filtered fragment library: Select fragments tailored for the recombination (remove pool X, deduplicate per subpocket, remove unfragmented ligands, remove all fragments that connect only to pool X, keep only fragment-like fragments, and filter for hinge-like AP fragments).
- `fragment_library_reduced/`: Reduced fragment library: Select a diverse set of fragments (per subpocket) for recombination starting from the filtered fragment library.
- `combinatorial_library/`: Combinatorial library based on the reduced fragment library.
- `external/`: Data from external resources.
