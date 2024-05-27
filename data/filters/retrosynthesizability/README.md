# retrosynthesizability data

- `retro.txt`: File containing the results that were already assembled in previous queries for fragment pairs or those that will be queried in future searches. It contains \[pair SMILES\]; \[child(ren) 1 SMILES\]; \[child(ren) 2 SMILES\]; \[plausibility/ies\] for every requested fragment pair.

**Note:** The ASKCOS results for the given KinFragLib data are already precomputed, thus for that ASKCOS does not need to be installed to successfully run this notebook. However, if the notebook is executed on new data, the ASKCOS needs to be installed beforehand. To install ASKCOS, please follow the installation given at [https://askcos-docs.mit.edu/](https://askcos-docs.mit.edu/guide/1-Introduction/1.1-Introduction.html). 