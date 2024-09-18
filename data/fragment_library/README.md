# KinFragLib: Full fragment library

The (full) fragment library resulting from the KinFragLib fragmentation procedure comprises about 3,000 fragments, 
which are the basis for exploring the subpocket-based chemical space of ligands co-crystallized with kinases. 

## Fragment library

Fragments are organized by the subpockets they occupy. Each fragment subpocket pool is stored in a SDF file:

	`AP.sdf`
	`FP.sdf`
	`GA.sdf`
	`SE.sdf`
	`B1.sdf`
	`B2.sdf`
	`X.sdf`

Each fragment contains the following information:

- 3D coordinates of the fragment's atoms as reported in the corresponding KLIFS complex structure file.
- `kinase`, `family`, and `group`: 
*Kinase* name, *family* and *group* of the kinase that the ligand (from which the fragment originates) was 
co-crystallized with.
- `complex_pdb`, `ligand_pdb`, `alt`, and `chain`: 
*PDB complex* and *ligand ID*, *alternate model* and *chain* for the KLIFS structure that the ligand 
(from which the fragment originates) was co-crystallized with.
- `atom.prop.subpocket`: Subpocket assignment for each of the fragment's atoms.
- `atom.prop.environment`: BRICS environment IDs for each of the fragment's atoms.

Please refer to `notebooks/kinfraglib/1_1_quick_start.ipynb` on how to load and work with this dataset.

## Original ligands

Original ligands that are composed of the fragments in the full fragment library are stored as CSV file:

    `original_ligands.csv`
    
Each ligand contains the following information:

- `kinase`, `family`, and `group`: 
*Kinase* name, *family* and *group* of the kinase that the ligand was co-crystallized with.
- `complex_pdb`, `ligand_pdb`, `alt`, and `chain`: 
*PDB complex* and *ligand ID*, *alternate model* and *chain* for the KLIFS structure that the ligand was co-crystallized with.
- `subpocket`: 
Subpockets that the ligand occupies.
- `ac_helix`:
aC-helix conformation for the KLIFS structure that the ligand was co-crystallized with.
- `smiles`:
Ligand's SMILES string.

Please refer to `notebooks/kinfraglib/2_1_fragment_analysis_original_ligands.ipynb` where this data is generated.
