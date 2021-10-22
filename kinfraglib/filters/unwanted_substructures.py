"""
Contains functions to filter out unwanted substructures provided by Brenk et al.
"""
import pandas as pd
from rdkit import Chem
from . import building_blocks


def get_brenk(fragment_library, DATA):
    """
    Getting the path to the unwanted substructures provided by Brenk et al. and filtering them out.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    DATA : str
        path to the csv file provided by Brenk

    Returns
    -------
    dict
        Containing
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional columns defining wether the fragment is accepted (1) or rejected (0).
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
        molecule = Chem.MolFromSmiles(row.smiles)
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
    d = dict()
    d["fragment_library"] = fragment_library_bool
    d["brenk"] = matches

    return d
