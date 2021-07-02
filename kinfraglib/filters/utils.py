"""
Helpful utility functions for custom_fragment_library.
"""
from rdkit import Chem
import pandas as pd


def save_smiles_wo_dummy(fragment_library, PATH_DATA):
    """
    Save smiles strings without dummy atoms for use in DataWarrior.

    Parameters
    ----------
    fragment_library : dict
        smiles series containing fragment smiles strings
    PATH_DATA : str
        Path to fragment library folder.
    -------

    """
    fragment_library = pd.concat(fragment_library).reset_index(drop=True)
    fragments = fragment_library["smiles"]
    path = str(str(PATH_DATA) + "/fragment_library/smiles_wo_dummy.sdf")
    writer = Chem.SDWriter(path)
    for fragment in fragments:
        writer.write(Chem.MolFromSmiles(fragment))
    writer.close()
