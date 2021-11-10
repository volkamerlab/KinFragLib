"""
Contains the SYBA filter.
"""
from syba.syba import SybaClassifier
from . import check
from . import utils
import pandas as pd


def calc_syba(fragment_library,
              cutoff=0,
              cutoff_criteria=">",
              column_name="bool_syba",
              query_type="mol",
              ):
    """
    Calculate the SYnthetic Bayesian Accessibility for each fragment and add a boolean column
    if the fragment is accepted for the defined cutoff or not

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets inculding all information
    cutoff : int
        defining the cutoff value for rejecting/accepting fragments. By default cutoff=0
    cutoff_criteria : str
        defining if the fragments values need to be >, <, >=, <=, == or != compared to the
        cutoff_value. By default cutoff_criteria=">"
    column_name : str
        defining the column name where the bool if the fragment is accepted (1) or rejected (0) is
        stored
    query_type : str
        "mol" or "smiles". Defining if the SYBA score gets predicted using the ROMol from the
        fragment library or the SMILES string. By default query_type = "mol".
    Returns
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional column (bool_syba) defining whether the fragment is accepted (1) or rejected (0)
        and the calculated SYBA score (syba) for each fragment.
    -------

    """
    sybas = []
    syba = SybaClassifier()
    syba.fitDefaultScore()
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    for subpocket in fragment_library.keys():
        pocketsyba = []
        fragment_library_df_subpocket = fragment_library_df.loc[
            fragment_library_df["subpocket"] == subpocket
        ]
        if query_type == "mol":
            for molecule in fragment_library_df_subpocket["ROMol"]:
                pocketsyba.append(syba.predict(mol=molecule))
            sybas.append(pocketsyba)
        elif query_type == "smiles":
            for smiles in fragment_library_df_subpocket["smiles"]:
                pocketsyba.append(syba.predict(smiles))
            sybas.append(pocketsyba)

    fragment_library_bool = check.accepted_rejected(
        fragment_library, sybas, cutoff, cutoff_criteria, column_name
    )
    fragment_library_bool = utils.add_values(fragment_library_bool, sybas, "syba")

    return fragment_library_bool
