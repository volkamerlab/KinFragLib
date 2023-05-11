# KinFragLib: Kinase-focused fragment library

[![GitHub Actions Build Status](https://github.com/volkamerlab/KinFragLib/workflows/CI/badge.svg)](https://github.com/volkamerlab/KinFragLib/actions?query=workflow%3ACI)

![KinFragLib workflow](./docs/img/toc_github_kinfraglib.png)

Please note that this repository is constantly updated. You can retrieve the repository state for the published KinFragLib paper in release [v1.0.0](https://github.com/volkamerlab/KinFragLib/releases/tag/v1.0.0).

## Table of contents

- [Description](#description)
- [Repository content](#repository-content)
- [Quick start](#quick-start)
- [Contact](#contact)
- [License](#license)
- [Citation](#citation)

## Repository content

This repository holds the following resources:  

1. Fragment library data and a link to the combinatorial library data.
2. *Quick start* notebook explaining how to load and use the fragment library.  
3. Notebooks 

    3.1. Notebooks covering the full analyses regarding the fragment and combinatorial libraries as described in 
          the corresponding paper.  
    3.2. WIP: Notebooks providing a custom filtering framework to reduce the fragment library size.
    
Please find detailed description of files in `data/` and `notebooks/` in the folders' `README` files.

## Description

**Exploring the kinase inhibitor space using subpocket-focused fragmentation and recombination**

Protein kinases play a crucial role in many cell signaling processes, 
making them one of the most important families of drug targets.
Fragment-based drug design has proven useful as one approach to develop novel kinase inhibitors. 
Usually, fragment-based methods follow a knowledge-driven approach, i.e., optimizing a focused set of fragments into 
molecular hits. 

We present here *KinFragLib*, a data-driven kinase-focused fragment library based on the structural kinome data 
retrieved from the [KLIFS](https://klifs.vu-compmedchem.nl) database.
Each kinase binding pocket (for DFG-in structures with non-covalent ligands) is automatically divided in *KinFragLib* 
into six subpockets, i.e. the adenine pocket (AP), front pocket (FP), solvent-exposed pocket (SE), gate area (GA) as 
well as back pocket 1 and 2 (B1 and B2), based on defined pocket-spanning residues.
Each co-crystallized ligand is fragmented using the BRICS algorithm and its fragments are assigned to the respective 
subpocket they occupy. 
Following this approach, a fragment library is created with respective subpocket pools. This fragment library enables 
an in-depth analysis of the chemical space of known kinase inhibitors, and can be used to enumerate recombined 
fragments in order to generate novel potential inhibitors.

WIP: Custom KinFragLib provides a pipeline to filter the fragments in KinFragLib checking for unwanted substructures (PAINS and Brenk et al.), lead-/drug-likeness (Rule of Three and QED), synthesizability (similarity to buyable building blocks and SYBA) and pairwise retrosynthesizability. Each filter can be (de-)activated and the parameters can be modified by the user to create a customized filtered fragment library. 

## Quick start

1. Clone this repository.

    ```bash
    git clone https://github.com/volkamerlab/KinFragLib.git
    ```

2. Create the `kinfraglib` conda environment. 

    ```bash
    # Change to KinFragLib directory
    cd /path/to/KinFragLib   
 
    # Create and activate environment
    conda env create -f environment.yml
    conda activate kinfraglib

    # Install the kinfraglib pip package
    cd ..
    pip install -e KinFragLib
    ```

3. Open the notebook `quick_start.ipynb` for an introduction on how to load and use the fragment library.

    ```bash
    # Change to KinFragLib directory (if you have not already)
    cd /path/to/KinFragLib
 
    # Start jupyter lab to explore the notebooks
    jupyter lab
    ```

## Contact

Please contact us if you have questions or suggestions.

* Open an issue on our GitHub repository: https://github.com/volkamerlab/KinFragLib/issues
* Or send us an email: andrea.volkamer@charite.de

We are looking forward to hearing from you!

## License

This resource is licensed under the [MIT](https://opensource.org/licenses/MIT) license, a permissive open source license.

## Citation

Sydow, D., Schmiel, P., Mortier, J., and Volkamer, A. KinFragLib: Exploring the Kinase Inhibitor Space Using Subpocket-Focused Fragmentation and Recombination. J. Chem. Inf. Model. 2020. https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00839

```bib
@article{doi:10.1021/acs.jcim.0c00839,
author = {Sydow, Dominique and Schmiel, Paula and Mortier, Jérémie and Volkamer, Andrea},
title = {KinFragLib: Exploring the Kinase Inhibitor Space Using Subpocket-Focused Fragmentation and Recombination},
journal = {Journal of Chemical Information and Modeling},
volume = {60},
number = {12},
pages = {6081-6094},
year = {2020},
doi = {10.1021/acs.jcim.0c00839},
note ={PMID: 33155465},
URL = {https://doi.org/10.1021/acs.jcim.0c00839}
}
```
