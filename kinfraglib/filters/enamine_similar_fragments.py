#!/usr/bin/env python

from __future__ import print_function
from rdkit import Chem
from rdkit import DataStructs
import numpy as np
from pathlib import Path
import argparse
from kinfraglib import utils, filters
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator


def read_enamine_sdf(path):
    """
    Read in Enamine Building Blocks from SDFile into a list.

    Parameters
    ----------
    path : str
        path to SDFile with all Enamine building blocks

    Returns
    -------
    obj
        list of RDKit molecules
    """

    mols = [mol for mol in Chem.SDMolSupplier(path) if if mol is not None]
    for mol in Chem.SDMolSupplier(path):
        if mol is not None:
            mols.append(mol)
    return mols


def calculate_fingerprints(mols):
    """Calculate RDKit fingerprints for given molecules

    Args:
        mols (list): List of RDKit molecules

    Returns:
        list: List of fingerprints 
    """
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fingerprints = [rdkit_gen.GetFingerprint(mol) for mol in mols]
    return fingerprints


def most_similar_fragment(fp_array, fp):
    """Calculate Tanimoto similarity and save most similar fragment

    Args:
        fp_array (list): List of fingerprints of Enamine building blocks
        fp (RDKit): Query fingerprint from KinFragLib

    Returns:
        float, int: Tanimoto similarity and index of most similar fragment
    """
    sim = 0
    ind = -1
    for i, x in enumerate(fp_array):
        s = DataStructs.TanimotoSimilarity(x, fp)
        if s > sim:
            sim = s
            ind = i
    return sim, ind


def find_most_similar_fragment(fragment_library, enamine_mols, file_path):
    """Find most similar Enamine fragment for each KinFragLib fragment not matching with Enamine

    Args:
        fragment_library (DataFrame): KinFragLib fragmentation library
        enamine_mols (list): list of RDKit molecules containing enamine building blocks
        file_path (str): path to output file
    """
    f = Chem.SDWriter(file_path)
    enamine_fps = calculate_fingerprints(enamine_mols)
    print(f"Calculated fingerprints")
    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    for frag in fragment_library_concat.ROMol:
        fp = calculate_fingerprints([frag])[0]
        sim, ind = most_similar_fragment(enamine_fps, fp)
        match = enamine_mols[ind]

        # rename enamine molecule to link to KinFragLib fragment  
        match.SetProp("_Name", str(frag.GetProp("_Name") + "_enamine"))
        # set Tanimoto similarity 
        match.SetProp("Similarity", str(sim))
        frag.SetProp("Similarity", str(sim))
        # write KinFragLib fragment to file 
        f.write(frag)
        # write matching Enamine fragment to file 
        f.write(match)
    f.close()


def apply_enamine_filter(fragment_library, path):
    """Apply Enamine building blocks filter

    Args:
        fragment_library (DataFrame): current fragment library per subpocket
        path (str): path to output file

    Returns:
        DataFrame: fragment library with building blocks filter added
    """
    fragment_library = filters.synthesizability.check_building_blocks(
        fragment_library, path
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
        help="output file name for most similar enamine sdf",
        required=True,
    )

    # parse cmd-line arguments
    args = parser.parse_args()

    PATH_ENAMINE = Path(args.enamine)
    PATH_FRAG_LIB = Path(args.fragmentlibrary)
    PATH_OUTPUT = Path(args.output)

    fragment_library = utils.read_fragment_library(PATH_FRAG_LIB)
    fragment_library = filters.prefilters.pre_filters(fragment_library)
    fragment_library = apply_enamine_filter(
        fragment_library, "data/filters/Enamine/Enamine_Building_Blocks.sdf",
    )
    print(f"Done reading in fragment library")
    # SDF contains all building blocks downloaded from enamine website
    enamine_mols = read_enamine_sdf(str(PATH_ENAMINE))
    print(f"Done reading enamine molecules")
    find_most_similar_fragment(
        fragment_library,
        enamine_mols,
        str(PATH_OUTPUT),
    )


if __name__ == "__main__":
    main()
