"""
Contains functions to apply the basic filtering steps
"""
import pandas as pd
from kinfraglib import utils

def pre_filters(fragment_library):
    """
    Gets a dict containing fragments organized in subpockets.
    Function
    - Removes pool X
    - Removes duplicates
    - Removes unfragmented ligands
    - Removes fragments only connecting to pool X
    And returns the fragment_library dict without those fragments
    
    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    
    Returns
    -------
    dict
        of pandas.DataFrames containing fragments organized in subpockets
    """
    #remove fragments in pool x
    [fragment_library.pop(x, None) for x in ['X']]
    
    #remove duplicates within subpockets
    fragment_library = _remove_duplicates(fragment_library)
    
    #remove fragments without dummy atoms (unfragmented ligands)
    fragment_library = _remove_unfragmented(fragment_library)
    
    #remove fragments only connecting to pool X
    fragment_library = _remove_connecting_only_x(fragment_library)
    
    #create a dictionary of fragments organie´zed in subpockets again
    fragment_library = _make_df_dict(fragment_library)

    return fragment_library

    
def _remove_duplicates(fragment_library):
    """
    removes duplicates from each subpocket of the fragment library
    
    Parameters
    ----------
    fragment_libray : list
        of fragments and their information
    
    Returns
    -------
    pandas DataFrame
        containing fragments without duplicates inside the subpockets
    """
    #remove duplicates
    fragment_library = pd.concat(fragment_library).reset_index(drop=True)
    fragment_library.groupby('subpocket', sort=False)
    # Get fragment count (by SMILES) per subpocket
    fragment_count = fragment_library.groupby(
    ['subpocket', 'smiles'],
    sort=False
    ).size()
    # Get first occurrence of SMILES per subpocket
    fragment_library = fragment_library.groupby(
    ['subpocket', 'smiles'],
    sort=False
    ).first()
    # Add fragment count to these representative fragments
    fragment_library['fragment_count'] = fragment_count
    fragment_library.reset_index(inplace=True)
    
    return fragment_library

def _remove_unfragmented(fragment_library):
    """
    removes fragments with no dummy atoms - unfragmented ligands.
    
    Parameters
    ----------
    fragment_libray : pandas DataFrame
        containing fragments
    
    Returns
    -------
    pandas DataFrame
        containing no unfragmented ligands
    """
    #remove fragments without dummy atoms (unfragmented)
    # Get fragments' (subpocket) connections 
    fragment_library['connections'] = utils.get_connections_by_fragment(fragment_library).connections
    # Unfragmented ligands?
    bool_unfragmented_ligands = fragment_library.connections.apply(lambda x: len(x) == 0)
    # Remove unfragmented ligands
    fragment_library = fragment_library[
    ~bool_unfragmented_ligands
    ].copy()
    
    return fragment_library

def _remove_connecting_only_x(fragment_library):
    """
    removes fragments which connect only to pool X
    
    Parameters
    ----------
    fragment_libray : pandas DataFrame
        containing fragments
    
    Returns
    -------
    pandas DataFrame
        containing without ligands that only connect to pool X
    """
    #remove fragments only connecting to pool x
    # Fragment connects only to pool X?
    bool_only_pool_x_connections = fragment_library.connections.apply(
    lambda x: all(  # All connections per fragment are X?
        [True if 'X' in i else False for i in x]  # Connections per fragment X or not?
    )
    )
    # Remove fragments that connect only to pool X
    fragment_library = fragment_library[
    ~bool_only_pool_x_connections
    ].copy()
    
    return fragment_library

def _make_df_dict(fragment_library):
    """
    Takes pandas DataFrame and creates dict to create the same format of the fragments as in the beginning.
    
    Parameters
    ----------
    fragment_libray : pandas DataFrame
        containing fragments
    
    Returns
    -------
    dict
        containing a pandas DataFrame for each subpocket
    """
    #reorder DataFrame into dict of pd.DataFrames again
    df = pd.DataFrame(fragment_library, columns=['ROMol', 'ROMol_dummy', 'ROMol_original', 'kinase', 'family', 'group',
       'complex_pdb', 'ligand_pdb', 'alt', 'chain', 'atom_subpockets',
       'atom_environments', 'smiles', 'smiles_dummy', 'subpocket'])
    fragment_library_dict = {}
    subpockets = fragment_library['subpocket'].unique()
    for subpocket in subpockets:
        fragment_library_dict[subpocket] = df[df.subpocket == subpocket]  
        fragment_library_dict[subpocket].reset_index(inplace = True, drop = True)
    fragment_library = fragment_library_dict
    return fragment_library