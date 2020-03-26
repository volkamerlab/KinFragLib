"""
Utility functions to work with the fragment library.
"""

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Lipinski
from rdkit.ML.Cluster import Butina
import seaborn as sns

SUBPOCKET_COLORS = {
    'AP': 'purple', 
    'FP': 'forestgreen', 
    'SE': 'c', 
    'GA': 'tab:orange', 
    'B1': 'tab:blue', 
    'B2': 'darkslateblue', 
    'X': 'grey'
}

def read_fragment_library(path_to_lib, remove_dummy=True):
    """
    Read fragment library from sdf files (one file per subpocket).
    
    Parameters
    ----------
    path_to_lib : str
        Path to fragment library folder.
    remove_dummy : bool
        Replace dummy atoms with hydrogens in fragments (default), or leave dummy atoms in fragments.
    
    
    Returns
    -------
    dict of pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, (values) for each subpocket (key).
    """
    
    # list of folders for each subpocket
    subpockets = ['AP', 'FP', 'SE', 'GA', 'B1', 'B2', 'X']
    
    data = {}

    # iterate over subpockets
    for subpocket in subpockets:

    	data[subpocket] = _read_subpocket_fragments(subpocket, path_to_lib, remove_dummy)
        
    return data

def _read_subpocket_fragments(subpocket, path_to_lib, remove_dummy=True):
    """
    Read fragments for input subpocket.
    
    Parameters
    ----------
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.
    path_to_lib : str
        Path to fragment library folder.
    remove_dummy : bool
        Replace dummy atoms with hydrogens in fragments (default), or leave dummy atoms in fragments.
    
    Returns
    -------
    pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, for input subpocket.
    """

    mol_supplier = Chem.SDMolSupplier(str(path_to_lib / f'{subpocket}.sdf'), removeHs=False)
        
    data = []

    for mol_raw in mol_supplier:
        
        if remove_dummy:
        
            # Replace dummy atoms with hydrogens in fragments
            dummy = Chem.MolFromSmiles('*')
            hydrogen = Chem.MolFromSmiles('[H]', sanitize=False)
            mol = AllChem.ReplaceSubstructs(mol_raw, dummy, hydrogen, replaceAll=True)[0]
            mol = Chem.RemoveHs(mol)  # Remove all hydrogens but explicit hydrogens
            
        else:
            
            mol = Chem.RemoveHs(mol_raw)
        
        # Generate SMILES for comparing fragments
        smiles = Chem.MolToSmiles(mol)
        
        # 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # kinase group
        data.append(
            [
                smiles,
                mol,
                mol.GetProp('kinase'),
                mol.GetProp('family'),
                mol.GetProp('group'),
                mol.GetProp('complex_pdb'),
                mol.GetProp('ligand_pdb'),
                mol.GetProp('alt'),
                mol.GetProp('chain'),
                mol.GetProp('atom.prop.subpocket'),
                mol.GetProp('atom.prop.environment')
            ]
        )

    return pd.DataFrame(
        data,
        columns='smiles fragment kinase family group complex_pdb ligand_pdb alt chain atom_subpockets atom_environments'.split()
    )

def most_common_fragments(fragments, top_x=50):
    """
    Get most common fragments.
    
    Parameters
    ----------
    fragments : pandas.DataFrame
        Fragment details, i.e. SMILES, kinase groups, and fragment RDKit molecules, for input subpocket.
        
    top_x : int
        Top x most common fragments.
        
    Returns
    -------
    tuple (list of rdkit.Chem.rdchem.Mol, pandas.Series)
        List of top x fragments (RDKit molecules) and frequence of top x fragments in subpocket (Series).
    """
    
    # Sort fragments by number of counts
    mols_count = fragments.smiles.value_counts()  # Sorted in descending order
    
    # Get RDKit Mol from SMILES
    mols = [Chem.MolFromSmiles(smiles) for smiles in mols_count.index]
    
    # N most common fragments
    return mols[:top_x], mols_count[:top_x]

def generate_fingerprints(mols):
    """
    Generate RDKit fingerprint from list of molecules.
    
    Parameters
    ----------
    mols : list of rdkit.Chem.rdchem.Mol
        List of molecules.
        
    Returns
    -------
    list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.
    """
    
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fingerprints = [rdkit_gen.GetFingerprint(mol) for mol in mols]
    
    return fingerprints

def cluster_molecules(fingerprints, cutoff=0.6):
    """
    Cluster fingerprints using the Butina algorithm.
    
    Parameters
    ----------
    fingerprints : list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.
    cutoff : float
        Distance cutoff Butina clustering.
        
    Returns
    -------
    list of tuple of int
        List of clusters, whereby each cluster is described by its cluster member IDs.
    """
    
    # Calculate Tanimoto distance matrix
    distance_matrix = _tanimoto_distance_matrix(fingerprints)
    
    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(
        distance_matrix,
        len(fingerprints),
        cutoff,
        isDistData=True
    )
    
    # Sort clusters by size
    clusters = sorted(clusters, key=len, reverse=True)
    
    # Get number of singleton clusters
    num_singletons = len([cluster for cluster in clusters if len(cluster) == 1])
    
    # Print details on clustering
    print("Number of fragments:", len(fingerprints))    
    print("Threshold: ", cutoff)
    print("Number of clusters: ", len(clusters))
    print("# clusters with only 1 compound: ", num_singletons)
    
    return clusters

def _tanimoto_distance_matrix(fingerprints):
    """
    Calculate distance matrix for list of fingerprints.
    
    Parameters
    ----------
    fingerprints : list of rdkit.DataStructs.cDataStructs.ExplicitBitVect
        List of fingerprints.
        
    Returns
    -------
    list of floats
        Distance matrix (a triangular similarity matrix in the form of a list)
    """
    
    fingerprints = list(fingerprints)
    distance_matrix = []
    
    for i in range(1,len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[:i])
        
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        distance_matrix.extend([1-x for x in similarities])
    
    return distance_matrix


def get_fragmented_ligand(fragment_library, complex_pdb, ligand_pdb):
    """
    Show fragments per subpocket for ligand by PDB ID.
    
    Parameters
    ----------
    fragment_library : dict of pandas.DataFrame
        Fragment details, i.e. SMILES, and fragment RDKit molecules, KLIFS and fragmentation details (values)
        for each subpocket (key).
    complex_pdb : str
        PDB ID for structure with ligand of interest.
    ligand_pdb : str
        PDB ID for ligand of interest.
    
    Returns
    -------
    PIL.PngImagePlugin.PngImageFile
        Fragmented ligand.
    """
    
    subpockets = ['SE', 'AP', 'GA', 'B1', 'B2', 'FP', 'X']  # order taken from paper Figure 4

    fragments = []

    for subpocket in subpockets:
        
        subpocket_fragments = fragment_library[subpocket]
        subpocket_fragments = subpocket_fragments[
            (subpocket_fragments.complex_pdb == complex_pdb) & (subpocket_fragments.ligand_pdb == ligand_pdb)
        ].copy()
        subpocket_fragments['subpocket'] = subpocket
        fragments.append(subpocket_fragments)

    fragmented_ligand = pd.concat(fragments)
    
    return fragmented_ligand


def draw_fragmented_ligand(fragment_library, complex_pdb, ligand_pdb, mols_per_row=6):
    """
    Show fragments per subpocket for ligand by PDB ID.
    
    Parameters
    ----------
    fragment_library : dict of pandas.DataFrame
        Fragment details, i.e. SMILES, and fragment RDKit molecules, KLIFS and fragmentation details (values)
        for each subpocket (key).
    complex_pdb : str
        PDB ID for structure with ligand of interest.
    ligand_pdb : str
        PDB ID for ligand of interest.
    
    Returns
    -------
    PIL.PngImagePlugin.PngImageFile
        Fragmented ligand.
    """
    
    fragmented_ligand = get_fragmented_ligand(fragment_library, complex_pdb, ligand_pdb)
    
    img = Draw.MolsToGridImage(
        fragmented_ligand.fragment.tolist(), 
        legends=fragmented_ligand.subpocket.tolist(), 
        molsPerRow=mols_per_row
    )
    
    return img


def _descriptors_from_mol(mol):
    """
    Get descriptors for a molecule, i.e. number of heavy atoms, logP, and number of hydrogen bond acceptors/donors.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    pd.Series
        Descriptors for input molecule.
    """

    size = mol.GetNumHeavyAtoms()
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    return pd.Series([size, logp, hbd, hba], index='size logp hbd hba'.split())


def descriptors_from_smiles(smiles):
    """
    Get descriptors for a set of SMILES.

    Parameters
    ----------
    smiles : pd.Series
        Set of SMILES.

    Returns
    -------
    pd.Series
        Descriptors for set of SMILES.
    """

    descriptors = pd.DataFrame(
        smiles.apply(
            lambda x: _descriptors_from_mol(Chem.MolFromSmiles(x))
        )
    )

    return descriptors


def descriptors_by_fragments(fragment_library):
    """
    Get physicochemical properties of fragment library, i.e. size (# heavy atoms), logP, hydrogen bond donors and acceptors.
    
    Parameters
    ----------
    fragment_library : dict of pandas.DataFrame
        Fragment details, i.e. SMILES, and fragment RDKit molecules, KLIFS and fragmentation details (values)
        for each subpocket (key).
        
    Returns
    -------
    pandas.DataFrame
        Properties of fragment library.
    """
    
    properties = {}

    for subpocket, fragments in fragment_library.items():
        properties[subpocket] = descriptors_from_smiles(fragments['smiles'].drop_duplicates())

    properties = pd.concat(properties).reset_index()

    properties.drop('level_1', axis=1, inplace=True)
    properties.rename(
        columns={
            'level_0': 'subpocket',
            'size': '# Heavy atoms',
            'logp': 'LogP',
            'hbd': '# HB donors',
            'hba': '# HB acceptors'
        },
        inplace=True
    )
    return properties


def _drug_likeness_from_mol(mol):
    """
    Get drug-likeness criteria for a molecule, i.e. molecular weight, logP, number of hydrogen bond acceptors/donors and
    accordance to Lipinski's rule of five.
    (Takes about 1s for 2000 mols.)

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule.

    Returns
    -------
    pd.Series
        Drug-likeness criteria for input molecule.
    """

    mw = 1 if Descriptors.ExactMolWt(mol) <= 500 else 0
    logp = 1 if Descriptors.MolLogP(mol) <= 5 else 0
    hbd = 1 if Lipinski.NumHDonors(mol) <= 5 else 0
    hba = 1 if Lipinski.NumHAcceptors(mol) <= 10 else 0
    lipinski = 1 if mw + logp + hbd + hba >= 3 else 0

    return pd.Series([mw, logp, hbd, hba, lipinski], index='mw logp hbd hba lipinski'.split())


def drug_likeness_from_smiles(smiles):
    """
    Get drug-likeness for a set of SMILES.

    Parameters
    ----------
    smiles : pd.Series
        Set of SMILES.

    Returns
    -------
    pd.Series
        Ratio of drug like molecules.
    """

    drug_likeness = pd.DataFrame(
        smiles.apply(
            lambda x: _drug_likeness_from_mol(Chem.MolFromSmiles(x))
        )
    )
    print(f'Number of molecules: {drug_likeness.shape[0]}')

    drug_likeness_ratio = round(drug_likeness.apply(sum) / len(drug_likeness) * 100, 0)

    return drug_likeness_ratio


def connections_by_fragment(fragment_library_concat):
    """
    For each fragment, extract connecting subpockets (e.g. ['FP', 'SE'] for subpocket 'AP') and define subpocket connections (e.g. ['AP=FP', 'AP=SE']). 
    
    Parameters
    ----------
    fragment_library_concat : pandas.DataFrame
        Fragment library data for one or mulitple subpockets.
        
    Returns
    -------
    pandas.DataFrame
        Fragment library data including connecting subpockets and connections.    
    """

    # For each fragment, extract connecting subpocket from atom_subpockets, e.g. ['FP', 'SE'] for subpocket 'AP'
    fragment_library_concat['connections'] = fragment_library_concat.apply(
        lambda x: [i for i in x.atom_subpockets.split() if i != x.subpocket], 
        axis=1
    )
    
    # Extract each connection (join connecting subpockets), e.g. ['AP=FP', 'AP=SE']
    fragment_library_concat['connections_name'] = fragment_library_concat.apply(lambda x: ["=".join(sorted([x.subpocket, i])) for i in x.connections], axis=1)

    return fragment_library_concat['kinase complex_pdb ligand_pdb atom_subpockets connections connections_name'.split()]


def connections_by_ligand(connections_by_fragment):
    """
    For each ligand, extract subpocket connections.
    
    Parameters
    ----------
    connections_by_fragment : pandas.DataFrame
        Fragment library data including connecting subpockets and connections (see connections_by_fragment() function).
        
    Returns
    -------
    pandas.DataFrame
        Ligands represented by fragment library with details on their subpocket connections. 
    """

    # Pool fragment connections by ligand
    connections_by_ligand = connections_by_fragment.groupby(['group', 'complex_pdb', 'ligand_pdb'])['connections_name'].sum()

    # Deduplicate connections (count each connection only once)
    connections_by_ligand = connections_by_ligand.apply(lambda x: set(x))

    return connections_by_ligand


def connections_count_by_ligand(connections_by_ligand):
    """
    Count subpocket connections (by type) across all ligands.
    
    Parameters
    ----------
    connections_by_ligand : pandas.DataFrame
        Ligands represented by fragment library with details on their subpocket connections (see connections_by_ligand() function). 
        
    Returns
    -------
    pandas.DataFrame
        Subpocket connections count and frequency across all ligands.
    """
    
    # For each ligand (row) count connection type (column)
    connection_matrix = pd.DataFrame({i: [] for i in connections_by_ligand.index}).transpose()

    for index, row in connections_by_ligand.iteritems():

        for connection in row:

            if connection not in connection_matrix.columns:
                connection_matrix[connection] = 0


            connection_matrix[connection][index] += 1
            
    # Count connection types per ligand
    connections_count = pd.DataFrame(
        {
            'count': connection_matrix.sum(), 
            'frequency': round(connection_matrix.sum() / connection_matrix.shape[0] * 100, 1)
        }
    ).sort_values('count', ascending=False)

    return connections_count


def fragment_similarity(fragment_library_concat, by='subpocket'):
    """
    Calculate similarities for all pairwise fragment combinations within each group.
    
    Parameters
    ----------
    fragment_library_concat : pandas.DataFrame
        Fragment library data for one or mulitple subpockets.
    by : str
        Group by which fragments are grouped and similarity is calculated, e.g. subpocket or kinase groups.
        
    Returns
    -------
    pandas.DataFrame
        Similarity values for each fragment pair per subpocket.
    """
    
    similarities_by_group = []

    for group_name, fragments in fragment_library_concat.groupby(by):

        smiles_deduplicated = fragments['smiles'].drop_duplicates()

        mols = smiles_deduplicated.apply(lambda x: Chem.MolFromSmiles(x))
        fingerprints = generate_fingerprints(mols)

        similarities = []

        for fp1, fp2 in combinations(fingerprints, 2):
            similarities.append(DataStructs.FingerprintSimilarity(fp1, fp2))
            
        similarities = pd.DataFrame(similarities)
        similarities.rename(columns={0: 'similarity'}, inplace=True)
        similarities[by] = group_name

        similarities_by_group.append(similarities)

    return pd.concat(similarities_by_group)


def plot_fragment_similarity(similarities_by_group):
    
    plt.figure(figsize=(8,8))
    
    try:
        ax = sns.boxplot(
            x=similarities_by_group.columns[1], 
            y=similarities_by_group.columns[0], 
            data=similarities_by_group, 
            palette=SUBPOCKET_COLORS
        )
    except KeyError:
        ax = sns.boxplot(
        x=similarities_by_group.columns[1], 
        y=similarities_by_group.columns[0], 
        data=similarities_by_group, 
        color='dodgerblue'
    )
    plt.ylabel('Tanimoto similarity', fontsize=18)
    plt.xlabel(similarities_by_group.columns[1], fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)