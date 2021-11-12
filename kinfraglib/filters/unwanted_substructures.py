"""
Contains functions to filter out unwanted substructures
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalogParams, FilterCatalog
from . import synthesizability


def get_pains(fragment_library):
    """
    Function to check fragments for PAINS structures.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information

    Returns
    -------
    fragment_library, matches: tuple(dict,dict)
        Containing
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional column (bool_pains) defining whether the fragment is accepted (1) or
            rejected (0).
            A pandas.DataFrame with the fragments and the names of the first PAINS structure found
            in the fragment.
    """
    # Code adapted from https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T003_compound_unwanted_substructures/talktorial.ipynb  # noqa: E501

    # initialize filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    # save fragment library as DataFrame
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    # search for PAINS
    matches = []
    clean = []
    accepted_bool = []
    for index, row in fragment_library_df.iterrows():
        molecule = Chem.MolFromSmiles(row.smiles)
        entry = catalog.GetFirstMatch(molecule)  # Get the first matching PAINS
        if entry is not None:
            # store PAINS information
            matches.append(
                {
                    "fragment": molecule,
                    "pains": entry.GetDescription().capitalize(),
                }
            )
            accepted_bool.append(0)
        else:
            # collect indices of molecules without PAINS
            clean.append(index)
            accepted_bool.append(1)
    # store fragment and pains structure found in the fragment
    matches = pd.DataFrame(matches)
    # add a boolean column if the fragment contains a pains structure
    fragment_library_bool = synthesizability._add_bool_column(
        fragment_library, accepted_bool, "bool_pains"
    )

    return fragment_library_bool, matches


def get_brenk(fragment_library, DATA):
    """
    Getting the path to the unwanted substructures provided by Brenk et al. and filtering them out.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information
    DATA : str
        path to the csv file provided by Brenk et al.

    Returns
    -------
    fragment_library, matches: tuple(dict,dict)
        Containing
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional column (bool_brenk) defining wether the fragment is accepted (1) or
            rejected (0).
            A pandas.DataFrame with the fragments, the substructures found and the substructure
            names
    """
    # Code adapted from https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T003_compound_unwanted_substructures/talktorial.ipynb # noqa: E501

    # read in csv file with unwanted substructure molecules
    substructures = pd.read_csv(DATA / "unwanted_substructures.csv", sep=r"\s+")
    substructures["rdkit_molecule"] = substructures.smarts.apply(Chem.MolFromSmarts)
    print(
        "Number of unwanted substructures in Brenk et al. collection:",
        len(substructures),
    )
    # save fragment library as DataFrame
    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)

    matches = []        # variable to store the matches (fragment and unwanted substructure found)
    clean = []          # variable to store the fragment indices without unwanted substructures
    rejected = []       # variable to store the fragment indices with unwanted substructures
    brenk_bool = []     # variable to store a bool for each fragment if unwanted substr. was found
    # iterate through rows of the fragment library Dataframe
    for index, row in fragment_library_df.iterrows():
        molecule = row.ROMol        # save molecule of fragment
        match = False
        # iterate through unwanted substructure molecules
        for _, substructure in substructures.iterrows():
            # check if the current fragment contains the unwanted substructure
            if molecule.HasSubstructMatch(substructure.rdkit_molecule):
                # if unwanted substructure is in fragment save fragment, unwanted substructure and
                # unwanted substructure name
                matches.append(
                    {
                        "fragment": molecule,
                        "substructure": substructure.rdkit_molecule,
                        "substructure_name": substructure["name"],
                    }
                )
                match = True        # set match to true
        if not match:       # fragment has no unwanted substructure
            clean.append(index)
            brenk_bool.append(1)
        else:       # unwanted substructure was found in fragment
            brenk_bool.append(0)
            rejected.append(index)
    # add unwanted substructures found to DataFrame
    matches = pd.DataFrame(matches)
    # add boolean column if an unwanted substructure was found to fragment library
    fragment_library_bool = synthesizability._add_bool_column(
        fragment_library, brenk_bool, "bool_brenk"
    )

    return fragment_library_bool, matches
