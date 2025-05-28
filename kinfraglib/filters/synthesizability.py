"""
Contains functions to filter for synthesizability
"""

import pandas as pd
from rdkit import Chem
from syba.syba import SybaClassifier
from . import check
from . import prefilters
from . import utils


def neutralize_atoms(mol):
    """Neutralize molecules (code taken from RDKit cookbook)

    Args:
        mol (RDKit Mol): input molecule

    Returns:
        RDKit mol: neutralized molecule
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def check_building_blocks(fragment_library, path_to_building_blocks):
    """
    Read in Enamine Building Blocks from SDFile created with filters/enamine_substructures.py
    and check if the fragment molecules are a substructure of building block molecules.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information
    path_to_building_blocks : str
        path to SDFile with overlapping building blocks is saved

    Returns
    -------
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional columns (bool_bb) defining whether the fragment is accepted (1), meaning found
        as a substructure in a building block, or rejected (0).
    """
    # save fragment library as DataFrame
    fragment_library_pre_filtered_df = pd.concat(fragment_library).reset_index(
        drop=True
    )
    bools_enamine = []
    # store Enamine Building Blocks from DatWarrioir SDF file
    bb_mols = _read_bb_sdf(path_to_building_blocks)
    print("Number of building blocks: %s" % len(bb_mols))

    # go through fragment library and Enamine Building Blocks and check if the fragments are
    # substructures of any Enamine Building Block loaded.
    for row in fragment_library_pre_filtered_df.itertuples():
        in_enamine = False
        frag_mol = row.ROMol
        for bb in bb_mols:
            if bb.HasSubstructMatch(frag_mol):
                in_enamine = True
                break
        if in_enamine:
            bools_enamine.append(1)
        else:
            bools_enamine.append(0)
    # add the boolean column if the fragment was found as a substrutcure of a Building Block
    fragment_library_bool = _add_bool_column(fragment_library, bools_enamine, "bool_bb")

    return fragment_library_bool


def _read_bb_sdf(path_to_building_blocks):
    """
    Read in Enamine Building Blocks from SDFile.

    Parameters
    ----------
    path_to_building_blocks : str
        path where SDFile with resulting Enamine building blocks is saved

    Returns
    -------
    list
        rdkit molecules of building blocks
    """
    enamine_bb = []
    # read in Enamine file with Enamine Building Blocks
    curpath = str(path_to_building_blocks)
    suppl = Chem.SDMolSupplier(curpath)
    # go through molecules from the read file and save it in a list
    for mol in suppl:
        enamine_bb.append(mol)
    return enamine_bb


def _add_bool_column(fragment_library, bool_list, column_name="bool"):
    """
    Adds a boolean column to the existing dict of pandas.DataFrames

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    bool_list : list
        containing boolean values
    column_name : str
        name the boolean column should be named

    Returns
    -------
    dict
        fragments organized in subpockets including boolean column
    """
    # save fragment library as a DataFrame
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    # add the boolean column to the DataFrame
    fragment_library_df[column_name] = pd.Series(
        bool_list, index=fragment_library_df.index
    )
    # create dict again with new column and return it
    fraglib = prefilters._make_df_dict(pd.DataFrame(fragment_library_df))
    return fraglib


def calc_syba(
    fragment_library,
    cutoff=0,
    cutoff_criteria=">",
    query_type="mol",
):
    """
    Calculate the SYnthetic Bayesian Accessibility (SYBA) for each fragment and add a boolean
    column if the fragment is accepted for the defined cutoff or not and a column with the
    calculated SYBA values.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information
    cutoff : int
        defining the cutoff value for rejecting/accepting fragments. By , cutoff=0
    cutoff_criteria : str
        defining if the fragment values need to be ">", "<", ">=", "<=", "==" or "!=" compared to
        the cutoff_value. By default, cutoff_criteria=">"
    query_type : str
        "mol" or "smiles". Defining if the SYBA score gets predicted using the ROMol from the
        fragment library or the SMILES string. By default, query_type = "mol".

    Returns
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional column (bool_syba) defining whether the fragment is accepted (1) or rejected (0)
        and the calculated SYBA score (syba) for each fragment.
    -------

    """
    sybas = []  # variable for storing the calculated SYBA values
    syba = SybaClassifier()  # loading the classifier to calculate the SYBA score
    syba.fitDefaultScore()  # fit the classifier to the default score
    # save fragment library as DataFrame
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    # iterate through subpockets
    for subpocket in fragment_library.keys():
        pocketsyba = []  # store syba values for every subpocket in a list
        # get all fragments from this subpocket
        fragment_library_df_subpocket = fragment_library_df.loc[
            fragment_library_df["subpocket"] == subpocket
        ]
        # calculate SYBA score for molecules if chosen
        if query_type == "mol":
            for molecule in fragment_library_df_subpocket["ROMol"]:
                pocketsyba.append(syba.predict(mol=molecule))
            sybas.append(pocketsyba)
        # calculate SYBA score for SMILES strings if chosen
        elif query_type == "smiles":
            for smiles in fragment_library_df_subpocket["smiles"]:
                pocketsyba.append(syba.predict(smiles))
            sybas.append(
                pocketsyba
            )  # add syba values from the subpocket to the syba list
    # add 'bool_syba' column to the fragment library
    fragment_library_bool = check.accepted_rejected(
        fragment_library, sybas, cutoff, cutoff_criteria, "bool_syba"
    )
    # add syba values to the fragment library
    fragment_library_bool = utils.add_values(fragment_library_bool, sybas, "syba")

    return fragment_library_bool
