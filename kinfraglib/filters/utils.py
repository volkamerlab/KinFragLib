"""
Helpful utility functions for custom_fragment_library.
"""
from rdkit import Chem
import pandas as pd


def save_fragments_wo_dummy(fragment_library, PATH_DATA):
    """
    Save fragments without dummy atoms in a .sdf file for use in DataWarrior.

    Parameters
    ----------
    fragment_library : dict
        fragment library organized in subpockets
    PATH_DATA : str
        Path where file should be saved.
    -------

    """
    fragment_library = pd.concat(fragment_library).reset_index(drop=True)
    fragments = fragment_library["ROMol"]
    path = str(str(PATH_DATA) + "/fragments_wo_dummy.sdf")
    writer = Chem.SDWriter(path)
    for fragment in fragments:
        writer.write(fragment)
    writer.close()


def add_values(fragment_library, values, colname):
    """
    Adding values to the fragment library.

    Parameters
    ----------
    fragment_library : dict
        fragment library organized in subpockets
    values : list
        containing the values that should be contained in the fragment library
    colname : str
        name of the new column with the values
    -------

    """
    pocket_num = 0
    for subpocket in fragment_library.keys():
        fragment_library[subpocket][colname] = values[pocket_num]
        pocket_num = pocket_num + 1
    return fragment_library
