# Kinase focused fragmentation library *KinFragLib*

## Exploring the chemical space of kinase inhibitors: Subpocket-based fragmentation for data-driven recombination

## Table of contents  

* [Introduction](#introduction)
* [Quick start](#quick-start)
* [Repository content](#repository-content)
* [Contact](#contact)
* [License](#license)
* [Citation](#citation)

### Introduction

(Back to [Table of contents](#table-of-contents).)

Protein kinases play a crucial role in many cell signaling processes, making them one of the most important families of drug targets.
Fragment-based drug design has proven useful as an approach to develop novel kinase inhibitors. However, fragment-based methods are usually limited to a knowledge-driven approach of optimizing a focused set of fragments. 
Here, we present a data-driven fragmentation and recombination approach instead. 
A novel computational fragmentation method was implemented, which splits known kinase inhibitors into fragments with respect to the subpockets that they occupy. Thereby, a fragment library with several pools, representing the subpockets, is created.
This fragment library enables an in-depth analysis of the chemical space of known kinase inhibitors, and is used to recombine fragments in order to generate novel potential inhibitors.

#### Fragmentation

For each input kinase-ligand complex, the kinase binding pocket is divided into six subpockets. The ligands are fragmented according to these subpockets, and a fragment library with several pools is created, where each pool corresponds to one subpocket and contains the fragments that were assigned to this subpocket.

<img src ="./docs/img/subpocket_centers.png" width = "600" align="left"> 
<br clear="all" />

#### Recombination

Every possible fragment recombination is enumerated in order to create a virtual combinatorial compound library. The fragments are reconnected only at the broken bonds, while preserving the original subpocket connection of each bond. 

### Quick start

(Back to [Table of contents](#table-of-contents).)

1. Clone this repository.

    ```bash
    git clone https://github.com/volkamerlab/KinFragLib.git
    ```

2. Create and activate `kffl` conda environment. 

    ```bash
    conda env create -f environment.yml
    conda activate kffl
    ```

3. Open `quick_start.ipynb` for an introduction on how to load and use the fragment library.

    ```bash
    cd /path/to/KinFragLib
    jupyter lab
    ```

### Repository content

(Back to [Table of contents](#table-of-contents).)

This repository holds (i) fragment library data and (ii) a *quick start* notebook explaining how to load and use the library alongside additional analysis notebooks.

    data/
    └── FragmentLibrary/
            ├── AP.sdf
            ├── FP.sdf
            ├── GA.sdf
            ├── SE.sdf
            ├── B1.sdf
            ├── B2.sdf
            └── X.sdf
    notebooks/
    ├── `quick_start.ipynb`
    ├── `cluster_most_common_fragments.ipynb`
    └── `TBA`
    

### Contact

(Back to [Table of contents](#table-of-contents).)

TBA

### License

(Back to [Table of contents](#table-of-contents).)

TBA

### Citation

(Back to [Table of contents](#table-of-contents).)

TBA
