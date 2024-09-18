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
    # save fragment library as a DataFrame
    fragment_library = pd.concat(fragment_library).reset_index(drop=True)
    fragments_mols = fragment_library["ROMol"]  # save molecules of the fragments
    path = str(str(PATH_DATA) + "/fragments_wo_dummy.sdf")  # path to save file
    # write molecules to file
    writer = Chem.SDWriter(path)
    for fragment_mol in fragments_mols:
        writer.write(fragment_mol)
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
    # iterate through subpockets
    pocket_num = (
        0  # helper variable to count which subpocket index is the current index
    )
    values = [val for val in values if len(val)] # colums that are not empty
    for subpocket in fragment_library.keys():
        # add value list with the current subpocket index to the fragment library
        fragment_library[subpocket][colname] = values[pocket_num]
        pocket_num = pocket_num + 1  # increase helper index for next subpocket
    return fragment_library
