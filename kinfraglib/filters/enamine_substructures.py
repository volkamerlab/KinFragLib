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
    # enamine_mols = []
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
                # enamine_mols.append(library.GetMol(indices[int(min_ind)]))
    f.close()
    # return enamine_mols


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


def main():
    
    parser = argparse.ArgumentParser()
    # add cmd-line arguments
    parser.add_argument('-e', '--enamine', type=str, help='file name of enamine building blocks sdf', required=True)
    parser.add_argument('-f', '--fragmentlibrary', type=str, help='path to fragment library', required=True)
    parser.add_argument('-o', '--output', type=str, help='output file name for matching building blocks sdf', required=True)
    
    # parse cmd-line arguments
    args = parser.parse_args()
    
    PATH_ENAMINE = Path(args.enamine)
    PATH_FRAG_LIB = Path(args.fragmentlibrary)
    PATH_OUTPUT = Path(args.output)

    fragment_library = utils.read_fragment_library(PATH_FRAG_LIB)
    fragment_library = filters.prefilters.pre_filters(fragment_library)
    print(f"Done reading in fragment library")
    # SDF contains all building blocks downloaded from enamine website
    enamine_library = read_enamine_sdf(
        PATH_ENAMINE
    )

    substructure_search(
        fragment_library,
        enamine_library,
        PATH_OUTPUT,
    )


if __name__ == "__main__":
    main()
