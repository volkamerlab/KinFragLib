"""
Contains functions to filter out unwanted substructures
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalogParams, FilterCatalog
from . import building_blocks


def get_pains(fragment_library):
    """
    Function to check fragments for PAINS structures.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets inculding all information

    Returns
    -------
    fragment_library, matches: tuple(dict,dict)
        Containing
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional column (bool_pains) defining wether the fragment is accepted (1) or
            rejected (0).
            A pandas.DataFrame with the fragments and the names of the first PAINS structure found
            in the fragment.
    """
    # Code adapted from https://github.com/volkamerlab/teachopencadd/blob/master/teachopencadd/talktorials/T003_compound_unwanted_substructures/talktorial.ipynb  # noqa: E501

    # initialize filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

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

    matches = pd.DataFrame(matches)

    fragment_library_bool = building_blocks._add_bool_column(
        fragment_library, accepted_bool, "bool_pains"
    )

    return fragment_library_bool, matches


def get_brenk(fragment_library, DATA):
    """
    Getting the path to the unwanted substructures provided by Brenk et al. and filtering them out.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets inculding all information
    DATA : str
        path to the csv file provided by Brenk

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

    substructures = pd.read_csv(DATA / "unwanted_substructures.csv", sep=r"\s+")
    substructures["rdkit_molecule"] = substructures.smarts.apply(Chem.MolFromSmarts)
    print(
        "Number of unwanted substructures in Brenk et al. collection:",
        len(substructures),
    )

    fragment_library_df = pd.concat(fragment_library).reset_index(drop=True)
    # search for PAINS
    matches = []
    clean = []
    rejected = []
    brenk_bool = []
    for index, row in fragment_library_df.iterrows():
        molecule = row.ROMol
        match = False
        for _, substructure in substructures.iterrows():
            if molecule.HasSubstructMatch(substructure.rdkit_molecule):
                matches.append(
                    {
                        "fragment": molecule,
                        "substructure": substructure.rdkit_molecule,
                        "substructure_name": substructure["name"],
                    }
                )
                match = True
        if not match:
            clean.append(index)
            brenk_bool.append(1)
        else:
            brenk_bool.append(0)
            rejected.append(index)

    matches = pd.DataFrame(matches)

    fragment_library_bool = building_blocks._add_bool_column(
        fragment_library, brenk_bool, "bool_brenk"
    )

    return fragment_library_bool, matches
