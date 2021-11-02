"""
Contains a function to filter fragments for PAINS structures.
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
    fragment_libray : dict
        fragments organized in subpockets inculding all information

    Returns
    -------
    dict
        Containing
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional columns defining wether the fragment is accepted (1) or rejected (0).
            A pandas.DataFrame with the fragments and the names of the PAINS structure found
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
    pains_dict = dict()
    pains_dict["fragment_library"] = fragment_library_bool
    pains_dict["pains"] = matches

    return pains_dict
