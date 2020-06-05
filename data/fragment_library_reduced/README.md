
# Reduced fragment library

The (full) fragment library resulting from the KinFragLib fragmentation procedure comprises of about 3,000 fragments, 
which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases
(see `data/fragment_library/`). 

In order to demonstrate how this library can be used for recombining ligands, we offer here a reduced fragment library 
based on the following filters:

1. Remove all fragments that are not useful in a recombination, i.e. duplicates, fragments in pool X, fragments without 
dummy atoms, and fragments with dummy atoms only connecting to pool X. 
Also remove all AP fragments that show no hydrogen bond donors and acceptors (not hinge-like).
2. Select a diverse set of fragments (per subpocket) for recombination to (i) save computational cost and 
(ii) avoid recombination of highly similar fragments.

Step 1 is necessary to focus on fragments meaningful for the recombination, 
whereas step 2 mainly aims to reduce computational costs during recombination.

Please refer to the notebook `notebooks/3_1_fragment_library_reduced.ipynb` to check how the data was generated and/or 
use it as a starting point for customized protocols.
