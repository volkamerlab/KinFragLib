"""
Contains the SYBA filter.
"""
from syba.syba import SybaClassifier

def calc_syba(smiles_series):
    """
    Go through SMILES Series ans calculate the SYnthetic Bayesian Accessibility.
    
    Parameters
    ----------
    smiles_series : pandas.Series
        smiles series containing fragment smiles strings
    
    
    Returns
    -------
    ??array??
        SYBA value for each SMILES
    """
    sybas = []
    syba = SybaClassifier()
    syba.fitDefaultScore()
    for smiles in smiles_series:
        sybas.append(syba.predict(smiles))
    return sybas