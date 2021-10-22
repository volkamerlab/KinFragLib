"""
Contains function to check for building blocks
"""
import pandas as pd
from . import prefilters
from rdkit import Chem


def check_building_blocks(fragment_library, path_to_building_blocks):
    """
    Read in Enamine Building Blocks from SDFile, compare fragments with building blocks.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    path_to_building_blocks : str
        path where SDFile with resulting building blocks from DataWarrior is saved

    Returns
    -------
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional columns defining wether the fragment is accepted (1) or rejected (0)
    """
    enamine_bb = _read_bb_sdf(path_to_building_blocks)
    fragment_library_pre_filtered_df = pd.concat(fragment_library).reset_index(
        drop=True
    )
    bools_enamine = []
    # store smiles as molecules from enamine
    bb_mols = []
    for buildingblock in enamine_bb[0]:
        bb_mols.append(Chem.MolFromSmiles(buildingblock))

    for row in fragment_library_pre_filtered_df.itertuples():
        if row.smiles in enamine_bb[0]:
            bools_enamine.append(1)
        else:
            in_enamine = False
            frag_mol = Chem.MolFromSmiles(row.smiles)
            for bb in bb_mols:
                if bb.HasSubstructMatch(frag_mol):
                    in_enamine = True
                    break
            if in_enamine:
                bools_enamine.append(1)
            else:
                bools_enamine.append(0)

    fragment_library_bool = _add_bool_column(fragment_library, bools_enamine, "bool_bb")

    return fragment_library_bool


def _read_bb_sdf(path_to_building_blocks):
    """
    Read in Enamine Building Blocks from SDFile.

    Parameters
    ----------
    path_to_building_blocks : str
        path where SDFile with resulting building blocks from DataWarrior is saved

    Returns
    -------
    ??
        smiles of building blocks
    """
    enamine_bb = []
    curpath = str(path_to_building_blocks)
    suppl = Chem.SDMolSupplier(curpath)
    curmols = []
    for mol in suppl:
        smiles = Chem.MolToSmiles(mol)
        curmols.append(smiles)
    enamine_bb.append(curmols)
    return enamine_bb


def _add_bool_column(fragment_library, bool_list, column_name="bool"):
    """
    Add boolean column to existing dict of pandas.DataFrames

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    bool_list : list
        containing boolean values
    column_name : str
        name the boolean column should be named

    Returns
    -------
    fragment_libray : dict
        fragments organized in subpockets inculding boolean column
    """
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    fragment_library_df[column_name] = pd.Series(
        bool_list, index=fragment_library_df.index
    )
    fraglib = prefilters._make_df_dict(pd.DataFrame(fragment_library_df))
    return fraglib
