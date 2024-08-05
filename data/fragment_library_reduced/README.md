# Reduced fragment library

The (full) fragment library resulting from the KinFragLib fragmentation procedure comprises 9505 fragments, which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases (see `data/fragment_library/`).

In order to demonstrate how this library can be used for recombining ligands, we offer here a reduced fragment library (727 fragments) based on the following filters:

1. Remove all fragments that are not useful in a recombination. Check `data/fragment_library_filtered/`.
2. Select a diverse set of fragments (per subpocket) for recombination to (i) save computational cost and (ii) avoid recombination of highly similar fragments.

Step 1 is necessary to focus on fragments tailored for the recombination, whereas step 2 mainly aims to reduce computational costs during recombination.

## Reduction steps

1. Perform Butina clustering per subpocket (`DISTANCE_CUTOFF`).
2. Select fragments from clusters: Most common fragment per cluster (if mulitple most common fragments, select the one most similar to cluster centroid)
  a. Select always one fragment per cluster, i.e. cluster centroid (`N_REPRESENTED_FRAGMENTS` = None).
  b. Select top X most common fragments per cluster, whereby X is depending on cluster size, e.g. 10 clustered fragments will be represented by one selected fragment (`N_REPRESENTED_FRAGMENTS` = 10).
3. Include fragments from singleton clusters? (`INCLUDE_SINGLETONS`)

## Selected reduction parameters

- `DISTANCE_CUTOFF` = 0.6
- `N_REPRESENTED_FRAGMENTS` = 10
- `INCLUDE_SINGLETONS` = True

Please refer to the notebook `notebooks/kinfraglib/3_1_fragment_library_reduced.ipynb` to check how the data was generated and/or use it as a starting point for customized protocols.