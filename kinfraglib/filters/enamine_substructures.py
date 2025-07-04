#!/usr/bin/env python

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
import numpy as np
from pathlib import Path
import argparse
from kinfraglib import utils, filters


def read_enamine_sdf(path):
    """
    Read in Enamine Building Blocks from SDFile into an RDKit substructure library.

    Parameters
    ----------
    path : str
        path to SDFile with all enamine building blocks

    Returns
    -------
    obj
        RDKit substructure library containing all molecules from given SDFile
    """
    # takes around 5min to load Enamine Building Blocks
    library = rdSubstructLibrary.SubstructLibrary(
        rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    )
    rdkit_errors = 0
    for mol in Chem.SDMolSupplier(path):
        try:
            library.AddMol(mol)
        except:
            rdkit_errors += 1
            continue

    print(f"{rdkit_errors} molecules could not be converted to RDKit molecule.")
    return library


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


def substructure_search(queries, library, path):
    """
    Performs a substructure match for each query molecule against a library of molecules and
    writes minimal enamine building block which matches a query to an SDF

    Parameters
    ----------
    queries : dict
        fragments organized in subpockets including all information
    library : RDKit.SubstructLibrary
        RDKit substructure library containing all Enamine building blocks
    path: str
        path to SDFile containing enamine building blocks which match at least one molecule from the query
    """
    f = Chem.SDWriter(path)
    for key in queries.keys():
        print(f"Substructure matching for {key} subpocket")
        for q in queries[key].ROMol:
            params = Chem.AdjustQueryParameters()
            params.adjustRingCount = True
            params.adjustRingChain = True
            # get all structure matches
            indices = library.GetMatches(
                Chem.AdjustQueryProperties(q, params), maxResults=-1
            )
            if len(indices) > 0:
                # get index to smallest fragment which matches
                min_ind = np.argmin(
                    [library.GetMol(ind).GetNumAtoms() for ind in indices]
                )
                f.write(library.GetMol(indices[int(min_ind)]))
    f.close()


def write_to_file(path, mols):
    """
    Write enamine building blocks which match at least one KinFragLib fragment.

    Parameters
    ----------
    path : str
        path to SDFile output file containing enamine building blocks
    mols : list
        list containing RDKit molecules of the enamine building blocks
    """
    with Chem.SDWriter(path) as w:
        for m in mols:
            w.write(m)


def apply_enamine_filter(fragment_library):
    fragment_library = filters.synthesizability.check_building_blocks(
        fragment_library, "data/filters/Enamine/enamine_substructures_neutralized.sdf"
    )
    for key in fragment_library.keys():
        fragment_library[key] = fragment_library[key].loc[
            fragment_library[key]["bool_bb"] == 0
        ]
    return fragment_library


def main():

    parser = argparse.ArgumentParser()
    # add cmd-line arguments
    parser.add_argument(
        "-e",
        "--enamine",
        type=str,
        help="file name of enamine building blocks sdf",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fragmentlibrary",
        type=str,
        help="path to fragment library",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="output file name for matching building blocks sdf",
        required=True,
    )

    # parse cmd-line arguments
    args = parser.parse_args()

    PATH_ENAMINE = Path(args.enamine)
    PATH_FRAG_LIB = Path(args.fragmentlibrary)
    PATH_OUTPUT = Path(args.output)

    fragment_library = utils.read_fragment_library(PATH_FRAG_LIB)
    fragment_library = filters.prefilters.pre_filters(fragment_library)
    print(f"Done reading in fragment library")
    # SDF contains all building blocks downloaded from enamine website
    enamine_library = read_enamine_sdf(str(PATH_ENAMINE))

    substructure_search(
        fragment_library,
        enamine_library,
        str(PATH_OUTPUT),
    )


if __name__ == "__main__":
    main()
