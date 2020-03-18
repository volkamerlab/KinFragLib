# Kinase focused fragmentation library *KinFragLib*

Subpocket-based fragmentation of kinase ligands from KLIFS

Protein kinases play a crucial role in many cell signaling processes, 
making them one of the most important families of drug targets.
Fragment-based drug design has proven useful as an approach to develop novel kinase inhibitors. 
Fragment-based methods usually follow the knowledge-driven approach of optimizing a focused set of fragments. 
We offer here a data-driven kinase-focused fragment library instead 
based on the structural kinome data from the [KLIFS](https://klifs.vu-compmedchem.nl) database.
Each kinase binding pocket (for DFG-in structures with non-covalent ligands) is divided into six subpockets, i.e.
the adenine pocket, front pocket, solvent-exposed pocket, gate area as well as back pocket 1 and 2.
The ligands are fragmented according to the occupied subpockets, 
creating a fragment library with respective subpocket pools. 
This fragment library enables an in-depth analysis of the chemical space of known kinase inhibitors, 
and is used to recombine fragments in order to generate novel potential inhibitors.


<img src ="./docs/img/subpocket_centers.png" width = "600" align="left"> 
<br clear="all" />

### Quick start

1. Clone this repository.

    ```bash
    git clone https://github.com/volkamerlab/KinFragLib.git
    ```

2. Create and activate the `kffl` conda environment. 

    ```bash
    conda env create -f environment.yml
    conda activate kffl
    ```

3. Open the notebook `quick_start.ipynb` for an introduction on how to load and use the fragment library.

    ```bash
    cd /path/to/KinFragLib
    jupyter lab
    ```

### Repository content

This repository holds (i) fragment library data and (ii) a *quick start* notebook explaining how to load and use the library alongside additional analysis notebooks.

    data/
    └── fragment_library/
            ├── AP.sdf
            ├── FP.sdf
            ├── GA.sdf
            ├── SE.sdf
            ├── B1.sdf
            ├── B2.sdf
            └── X.sdf
    notebooks/
    ├── quick_start.ipynb
    ├── cluster_most_common_fragments.ipynb
    └── TBA
    

### Contact

Please contact us if you have questions or suggestions.

* Open an issue on our GitHub repository: https://github.com/volkamerlab/KinFragLib/issues
* Or send us an email: andrea.volkamer@charite.de

We are looking forward to hearing from you!

### License

TBA

### Citation

TBA
