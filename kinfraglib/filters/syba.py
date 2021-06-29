"""
Contains the SYBA filter.
"""
from syba.syba import SybaClassifier
from . import check
import pandas as pd

def calc_syba(fragment_library, cutoff=0, cutoff_criteria = "<", column_name = "bool_syba"):
    """
    Go through SMILES Series and calculate the SYnthetic Bayesian Accessibility.
    
    Parameters
    ----------
    fragment_library : dict
        smiles series containing fragment smiles strings
    cutoff : int
        defining the cutoff value for rejecting/accepting fragments
    cutoff_criteria : str
        defining if the fragments values need to be >, <, >=, <=, == or != compared to the cutoff_value
    column_name : str
        defining the column name where the bool if the fragment is accepted (1) or rejected (0) is stored
    Returns
    dict
        Containing
            A pandas.DataFrame with accepted fragments and their information.
            A pandas.DataFrame with rejected fragments and their information
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an additional columns defining
            wether the fragment is accepted (1) or rejected (0)
    -------
    
    """
    sybas = []
    syba = SybaClassifier()
    syba.fitDefaultScore()
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    for subpocket in fragment_library.keys():
        pocketsyba = []
        fragment_library_df_subpocket = fragment_library_df.loc[fragment_library_df['subpocket'] == subpocket]
        for smiles in fragment_library_df_subpocket['smiles']:
            pocketsyba.append(syba.predict(smiles))
        sybas.append(pocketsyba)
    
    df_accepted, df_rejected, fragment_library_bool = check.accepted_rejected(fragment_library, sybas, cutoff, cutoff_criteria, column_name)
    d = dict()
    d['df_accepted'] = df_accepted
    d['df_rejected'] = df_rejected
    d['fragment_library'] = fragment_library_bool
    d['sybas'] = sybas
    return d