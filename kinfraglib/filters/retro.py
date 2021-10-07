from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Chem import AllChem
from functools import reduce
from . import brics_rules
import multiprocessing as mp
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import requests
import copy


def pairwise_retrosynthesis(fragment_library_filtered):
    res = get_valid_pairs(fragment_library_filtered)
    valids = checkvalid(res, fragment_library_filtered)
    bonds = get_bonds(valids, res, fragment_library_filtered)
    pair_df = get_pairs(valids, bonds, fragment_library_filtered)
    # only for testing with subset
    pair_df = pair_df[0:10]
    num_cpu = mp.cpu_count()
    # create list of smiles as parallel computing cannot handle molecules
    pairs_smiles = []
    for mol in pair_df["pair"]:
        pairs_smiles.append(Chem.MolToSmiles(mol))
    df_split = np.array_split(pairs_smiles, num_cpu)
    para_res = Parallel(n_jobs=num_cpu)(
        delayed(call_retro_parallel)(split) for split in df_split
    )
    para_result = pd.concat(para_res)
    children_list = []
    for i, row in para_result.iterrows():
        for num_children in range(len(row["child 1"])):
            children_list.append(row["child 1"][num_children])
            children_list.append(row["child 2"][num_children])
    children_list = set(children_list)

    children_mols = []
    children_smiles = []
    for smile in children_list:
        if smile is not None:
            children_mols.append(Chem.MolFromSmiles(smile))
            children_smiles.append(smile)
        else:
            children_mols.append(None)
            children_smiles.append(None)

    pairs_frags_smiles = []
    frag1 = []
    frag2 = []
    pair = []
    for fragids in pair_df["fragment ids"]:
        frag1.append(
            Chem.MolToSmiles(
                fragment_library_filtered[fragids[0].split("_")[0]]["ROMol"][
                    int(fragids[0].split("_")[1])
                ]
            )
        )
        frag2.append(
            Chem.MolToSmiles(
                fragment_library_filtered[fragids[1].split("_")[0]]["ROMol"][
                    int(fragids[1].split("_")[1])
                ]
            )
        )
    for pairmol in pair_df["pair"]:
        pair.append(Chem.MolToSmiles(pairmol))
    pairs_frags_smiles = pd.DataFrame(
        list(zip(pair_df["fragment ids"], frag1, frag2, pair)),
        columns=("fragment ids", "fragment 1", "fragment 2", "pair"),
    )
    pairs_frags_smiles

    res_df = compare_mols(para_result, pairs_frags_smiles)
    countfrag, fraglib_filtered = retro_fragments(res_df, fragment_library_filtered)

    # Todo: write function to make res_df a pandas df with molecules
    mol_df = get_mol_df(res_df)
    return fraglib_filtered, mol_df, countfrag


def get_mol_df(res_df):
    frag1_mol = []
    frag2_mol = []
    pair_mol = []
    child1_mol = []
    child2_mol = []
    for i, row in res_df.iterrows():
        frag1_mol.append(Chem.MolFromSmiles(row["fragment 1"]))
        frag2_mol.append(Chem.MolFromSmiles(row["fragment 2"]))
        pair_mol.append(Chem.MolFromSmiles(row["pair"]))
        if row["child 1"] is not None:
            child1_mol.append(Chem.MolFromSmiles(row["child 1"]))
            child2_mol.append(Chem.MolFromSmiles(row["child 2"]))
        else:
            child1_mol.append(None)
            child2_mol.append(None)

    mol_df = pd.DataFrame(
        list(
            zip(
                res_df["fragment ids"],
                frag1_mol,
                frag2_mol,
                pair_mol,
                child1_mol,
                child2_mol,
                res_df["plausibility"],
            )
        ),
        columns=(
            "fragment ids",
            "fragment 1",
            "fragment 2",
            "pair",
            "child 1",
            "child 2",
            "plausibility",
        ),
    )

    return mol_df


def compare_mols(para_result, pairs_frags_smiles):
    para_res = para_result.copy(deep=True)
    para_res.set_index("pair", inplace=True)
    smiles_list = list(pairs_frags_smiles["fragment 1"])
    smiles_list.extend(list(pairs_frags_smiles["fragment 2"]))
    children_list = []
    for i, row in para_res.iterrows():
        for num_children in range(len(row["child 1"])):
            children_list.append(row["child 1"][num_children])
            children_list.append(row["child 2"][num_children])
    smiles_list.extend(list(children_list))
    smiles_list = set(smiles_list)  # unique list of all smiles (fragments and children)
    mols = get_mol(smiles_list)
    mols.set_index(
        "smiles", inplace=True
    )  # get mol from specific smiles mols.loc['Cc1cc(N)[nH]n1']['mol']
    # dataframe for result which frags are matching
    column_names = [
        "fragment ids",
        "fragment 1",
        "fragment 2",
        "pair",
        "child 1",
        "child 2",
        "plausibility",
    ]
    result_df = pd.DataFrame(columns=column_names)

    for i, row in pairs_frags_smiles.iterrows():
        cur_pair_smiles = row["pair"]
        cur_frag1_smiles = row["fragment 1"]
        frag1_mol = mols.loc[cur_frag1_smiles]["mol"]
        cur_frag2_smiles = row["fragment 2"]
        frag2_mol = mols.loc[cur_frag2_smiles]["mol"]
        frag_ids = row["fragment ids"]
        cur_children1_smiles = para_res.loc[cur_pair_smiles]["child 1"]
        cur_children2_smiles = para_res.loc[cur_pair_smiles]["child 2"]
        cur_probs = para_res.loc[cur_pair_smiles]["plausibility"]
        # go through children lists and compare
        for num_cur_smiles in range(len(cur_children1_smiles)):
            child1_smiles = cur_children1_smiles[num_cur_smiles]
            child2_smiles = cur_children2_smiles[num_cur_smiles]
            child1_mol = mols.loc[child1_smiles]["mol"]
            child2_mol = mols.loc[child2_smiles]["mol"]
            if child1_mol is not None and child2_mol is not None:
                if child1_mol.HasSubstructMatch(
                    frag1_mol
                ) and child2_mol.HasSubstructMatch(frag2_mol):
                    result_df = result_df.append(
                        {
                            "fragment ids": frag_ids,
                            "fragment 1": cur_frag1_smiles,
                            "fragment 2": cur_frag2_smiles,
                            "pair": cur_pair_smiles,
                            "child 1": child1_smiles,
                            "child 2": child2_smiles,
                            "plausibility": cur_probs[num_cur_smiles],
                        },
                        ignore_index=True,
                    )
                elif child1_mol.HasSubstructMatch(
                    frag2_mol
                ) and child2_mol.HasSubstructMatch(frag1_mol):
                    result_df = result_df.append(
                        {
                            "fragment ids": frag_ids,
                            "fragment 1": cur_frag1_smiles,
                            "fragment 2": cur_frag2_smiles,
                            "pair": cur_pair_smiles,
                            "child 1": child2_smiles,
                            "child 2": child1_smiles,
                            "plausibility": cur_probs[num_cur_smiles],
                        },
                        ignore_index=True,
                    )
            else:
                result_df = result_df.append(
                    {
                        "fragment ids": frag_ids,
                        "fragment 1": cur_frag1_smiles,
                        "fragment 2": cur_frag2_smiles,
                        "pair": cur_pair_smiles,
                        "child 1": None,
                        "child 2": None,
                        "plausibility": 0,
                    },
                    ignore_index=True,
                )

    return result_df


def retro_fragments(retro_df, fragment_library):
    fraglib = copy.deepcopy(fragment_library)
    # get list of fragment ids
    all_frags = []
    frag_ids = []
    for i, row in retro_df.iterrows():
        if row["plausibility"] != 0:
            frag_ids.append(retro_df["fragment ids"][i][0])
            frag_ids.append(retro_df["fragment ids"][i][1])
    all_frags = pd.DataFrame(frag_ids, columns=["ids"])
    # count number of frags
    counts = all_frags.groupby("ids").size()

    # go through all subpockets and fragments and add number of
    # contributions to retrosynth. pathways
    for subpocket in fraglib.keys():
        count_frags = []
        for i in range(0, len(fraglib[subpocket])):
            if hasattr(counts, str(subpocket + "_" + str(i))):
                attribute = str(subpocket + "_" + str(i))
                num_counts = getattr(counts, attribute)
                count_frags.append(num_counts)

            else:
                count_frags.append(0)
        fraglib[subpocket]["retro_count"] = count_frags

    return counts, fraglib


def get_mol(smiles_list):
    mols = []
    smiles = []
    for smile in smiles_list:
        if smile is not None:
            mols.append(Chem.MolFromSmiles(smile))
            smiles.append(smile)
        else:
            mols.append(None)
            smiles.append(None)
    df = pd.DataFrame(list(zip(mols, smiles)), columns=("mol", "smiles"))
    return df


def call_retro_parallel(pair_smiles):
    """
    One step retrosynthesis using ASKCOS for all valid build pairs of fragments. 
    Saving the plausibility and the children that can build this pair according to retrosynthetic
    analysis.

    Parameters
    ----------
    pair_smiles : numpy array
        containing SMILES strings of pairs build by fragments

    Returns
    -------
    pandas DataFrame
        containing the pair, the children building this pair and their plausibility

    """
    pairs = []
    children1 = []
    children2 = []
    plausibilities = []
    for smile in pair_smiles:
        pairs.append(smile)
        cur_children1 = []
        cur_children2 = []
        cur_plausibilities = []
        HOST = "https://askcos.mit.edu/"
        params = {
            "smiles": smile,  # required
            # optional with defaults shown
            "max_depth": 1,  # maximum number of reaction steps
            "max_branching": 25,  # ?max number of branches are looked at to find "best"?
            "expansion_time": 20,  # how long the expansion can run
            "max_ppg": 100,  # maximum price per gram
            "template_count": 100,
            # "max_cum_prob"
            # which common probability reached until no more templates are used
            "max_cum_prob": 0.995,
            # "chemical_property_logic"
            # molecules are buyable or not, can be 'none' (only price relevant),
            # 'and' (price and heavy atoms constraint) or
            # 'or' (one of both constraints is relevant)
            "chemical_property_logic": "none",
            # max heavy atom contraints if 'and' or 'or' is used in 'chemical_property_logic'
            "max_chemprop_c": 0,
            "max_chemprop_n": 0,
            "max_chemprop_o": 0,
            "max_chemprop_h": 0,
            # want to use popular chemicals as reasonable stopping points?
            "chemical_popularity_logic": "none",
            "min_chempop_reactants": 5,  # min frequence as popular reactant
            "min_chempop_products": 5,  # min frequence as popular prouct
            "filter_threshold": 0.75,
            "return_first": "true",  # default is false
        }
        resp = requests.get(HOST + "/api/treebuilder/", params=params, verify=False)
        # res.append(resp.json())
        retro = resp.json()

        if (len(retro["trees"])) > 0:
            for num_tree in range(0, len(retro["trees"])):
                if len(retro["trees"][num_tree]["children"][0]["children"]) == 2:
                    plausibility = retro["trees"][0]["children"][0]["plausibility"]
                    child1 = retro["trees"][num_tree]["children"][0]["children"][0]["smiles"]
                    child2 = retro["trees"][num_tree]["children"][0]["children"][1]["smiles"]
                    cur_children1.append(child1)
                    cur_children2.append(child2)
                    cur_plausibilities.append(plausibility)

        else:
            cur_children1.append(None)
            cur_children2.append(None)
            cur_plausibilities.append(0)
        children1.append(cur_children1)
        children2.append(cur_children2)
        plausibilities.append(cur_plausibilities)
    res = pd.DataFrame(
        list(zip(pairs, children1, children2, plausibilities)),
        columns=["pair", "child 1", "child 2", "plausibility"],
    )
    return res


def construct_ligand(fragment_ids, bond_ids, fragment_library):
    """
    *copied and adapted from kinase_focused_fragment_library*
    Construct a ligand by connecting multiple fragments based on a Combination object
    Parameters
    ----------
    fragment_ids: list of str
        Fragment IDs of recombined ligand, e.g. `["SE_2", "AP_0", "FP_2"]` (`<subpocket>_<fragment
        index in subpocket pool>`).
    bond_ids : list of list of str
        Bond IDs of recombined ligand, e.g. `[["FP_6", "AP_10"], ["AP_11", "SE_13"]]`:
        Atom (`<subpocket>_<atom ID>`) pairs per fragment bond.
    fragment_library : dict of pandas.DataFrame
        SMILES and RDKit molecules for fragments (values) per subpocket (key).
    Returns
    -------
    ligand: rdkit.Chem.rdchem.Mol or None
        Recombined ligand (or None if the ligand could not be constructed)
    """

    fragments = []
    for fragment_id in fragment_ids:

        # Get subpocket and fragment index in subpocket
        subpocket = fragment_id.split("_")[0]
        fragment_index = int(fragment_id.split("_")[1])
        fragment = fragment_library[subpocket].ROMol_original[fragment_index]

        # Store unique atom identifiers in original molecule (important for recombined ligand
        # construction based on atom IDs)
        fragment = Chem.RemoveHs(fragment)
        for i, atom in enumerate(fragment.GetAtoms()):
            fragment_atom_id = f"{subpocket}_{i}"
            atom.SetProp("fragment_atom_id", fragment_atom_id)
            atom.SetProp("fragment_id", fragment.GetProp("complex_pdb"))
        fragment = PropertyMol(fragment)

        # Append fragment to list of fragments
        fragments.append(fragment)

    # Combine fragments using map-reduce model
    combo = reduce(Chem.CombineMols, fragments)

    bonds_matching = True
    ed_combo = Chem.EditableMol(combo)
    replaced_dummies = []

    # for bond in bond_ids:

    dummy_1 = next(
        atom
        for atom in combo.GetAtoms()
        if atom.GetProp("fragment_atom_id") == bond_ids[0]
    )
    dummy_2 = next(
        atom
        for atom in combo.GetAtoms()
        if atom.GetProp("fragment_atom_id") == bond_ids[1]
    )
    atom_1 = dummy_1.GetNeighbors()[0]
    atom_2 = dummy_2.GetNeighbors()[0]

    # check bond types
    bond_type_1 = combo.GetBondBetweenAtoms(
        dummy_1.GetIdx(), atom_1.GetIdx()
    ).GetBondType()
    bond_type_2 = combo.GetBondBetweenAtoms(
        dummy_2.GetIdx(), atom_2.GetIdx()
    ).GetBondType()
    if bond_type_1 != bond_type_2:
        bonds_matching = False
        print("Bonds not matching")

    ed_combo.AddBond(atom_1.GetIdx(), atom_2.GetIdx(), order=bond_type_1)

    replaced_dummies.extend([dummy_1.GetIdx(), dummy_2.GetIdx()])

    # Do not construct this ligand if bond types are not matching
    if not bonds_matching:
        return

    # Remove replaced dummy atoms
    replaced_dummies.sort(reverse=True)
    for dummy in replaced_dummies:
        ed_combo.RemoveAtom(dummy)

    ligand = ed_combo.GetMol()

    # Replace remaining dummy atoms with hydrogens
    du = Chem.MolFromSmiles("*")
    h = Chem.MolFromSmiles("[H]", sanitize=False)
    ligand = AllChem.ReplaceSubstructs(ligand, du, h, replaceAll=True)[0]
    try:
        ligand = Chem.RemoveHs(ligand)
    except ValueError:
        print(Chem.MolToSmiles(ligand))
        return

    # Clear properties
    for prop in ligand.GetPropNames():
        ligand.ClearProp(prop)
    for atom in ligand.GetAtoms():
        atom.ClearProp("fragment_atom_id")

    # Generate 2D coordinates
    AllChem.Compute2DCoords(ligand)

    return ligand


def get_bonds(valids, data, fragment_library):
    """
    Function for getting the ccorresponding bond type to the connections of fragment pairs.

    Parameters
    ----------
    valids : list
        list of lists containing fragment id pairs of mathcing pairs
    
    data : dict
        fragment library prepared for building valid pairs

    fragment_libray : dict
        fragments organized in subpockets inculding all information


    Returns
    -------
    list
        list of lists containing fragment ids of pairs and corresponding bond type

    """

    bonds = (
        []
    )  # store bonds of valid matching pairs as atom IDs where connection is formed
    # go through all valid pairs
    for valid in valids:
        bond = []
        for val in valid:
            # load fragments that should get connected
            subpocket1 = val[0].split("_")[0]
            fragment1_index = int(val[0].split("_")[1])
            fragment1 = fragment_library[subpocket1]["ROMol_original"][fragment1_index]
            # remove Hs before finding bonds otherwise bond ids not correct because for combining
            # molecules without Hs are used
            fragment1 = Chem.RemoveHs(fragment1)

            subpocket2 = val[1].split("_")[0]
            fragment2_index = int(val[1].split("_")[1])
            fragment2 = fragment_library[subpocket2]["ROMol_original"][fragment2_index]
            # remove Hs before finding bonds
            fragment2 = Chem.RemoveHs(fragment2)

            # i = 0
            bond1_id = None
            bond2_id = None

            data1 = data[subpocket1][fragment1_index]
            # get corresponding connection to load environment, bond type and neighboring subpocket
            for i in range(0, len(data1.ports)):
                environment1 = data1.ports[i].environment
                bond_type1 = data1.ports[i].bond_type
                neighbor1 = data1.ports[i].neighboring_subpocket

                data2 = data[subpocket2][
                    fragment2_index
                ]  # for matching fragment also get the connection data
                for j in range(0, len(data2.ports)):
                    environment2 = data2.ports[j].environment
                    bond_type2 = data2.ports[j].bond_type
                    neighbor2 = data2.ports[j].neighboring_subpocket

                    # check again if BRICS bond, bond types and subpockets are matching for
                    #  a connection
                    if (
                        brics_rules.is_brics_bond(environment1, environment2)
                        and bond_type1 == bond_type2
                        and subpocket2 == neighbor1
                        and subpocket1 == neighbor2
                    ):
                        # get atom indices where connection is build
                        for atom in fragment1.GetAtoms():
                            atom_symbol = atom.GetSymbol()
                            if atom_symbol == "*":
                                bond1_id = subpocket1 + "_" + str(atom.GetIdx())

                        for atom2 in fragment2.GetAtoms():
                            atom_symbol2 = atom2.GetSymbol()
                            if atom_symbol2 == "*":
                                bond2_id = subpocket2 + "_" + str(atom2.GetIdx())

            bond.append(
                [bond1_id, bond2_id, bond_type1]
            )  # save atom indices and bond type for building the connection
        bonds.append(bond)
    return bonds


def get_pairs(valids, bonds, fragment_library):
    """
    Function to get built pairs from fragments, corresponding fragments and fragment ids.

    Parameters
    ----------
    valids : list
        list of lists containing fragment id pairs of mathcing pairs

    bonds : list
        list of lists containing fragment id pairs and corresponding bond type

    fragment_libray : dict
        fragments organized in subpockets inculding all information

    Returns
    -------
    pandas DataFrame
        containing fragment ids, fragments building pairs and paired molecules

    """
    pairs = []
    frags1 = []
    frags2 = []
    ids = []
    for i in range(0, len(valids)):
        for j in range(0, len(valids[i])):
            frag1 = fragment_library[valids[i][j][0].split("_")[0]][
                "ROMol_dummy"
            ][int(valids[i][j][0].split("_")[1])]
            frag2 = fragment_library[valids[i][j][1].split("_")[0]][
                "ROMol_dummy"
            ][int(valids[i][j][1].split("_")[1])]

            frags1.append(frag1)
            frags2.append(frag2)

            pair = construct_ligand(
                valids[i][j], bonds[i][j], fragment_library
            )
            pairs.append(pair)
            ids.append(valids[i][j])

    return pd.DataFrame(
        {"fragment ids": ids, "fragment1": frags1, "fragment2": frags2, "pair": pairs}
    )


def checkvalid(data, fragment_library):

    """
    Function for checking if the fragment pairs are valid and can build connections.

    Parameters
    ----------
    data : dict
         fragment library prepared for building valid pairs

    fragment_libray : dict
        fragments organized in subpockets inculding all information

    Returns
    -------
    list
        list of lists containing fragment id pairs of mathcing pairs

    """

    matches = []  # save matching fragment pairs
    # iterate through subpockets
    for subpocket in fragment_library.keys():
        # iterate through fragments in subpockets
        for fragment in data[subpocket]:
            fragment_id1 = (
                fragment.frag_id
            )  # store fragment ID of first fragment in pair
            # go through atom connnections and check neighbors, bond type and environment
            for i in range(0, len(fragment.ports)):
                neighbor = fragment.ports[i].neighboring_subpocket
                bond_type = fragment.ports[i].bond_type
                environment = fragment.ports[i].environment
                match = []  # store current matching fragment pair
                for frag2 in data[neighbor]:
                    fragment_id2 = frag2.frag_id  # store fragment ID of second fragment
                    for i in range(0, len(frag2.ports)):
                        # check environment type, subpocket, bond type
                        environment_match = brics_rules.is_brics_bond(
                            environment, frag2.ports[i].environment
                        )  # check if BRICS environments are able to form connection
                        # if subpocket is adjacent, bond type is eqal and environments are
                        # matching, add as valid matching pair
                        if (
                            frag2.ports[i].neighboring_subpocket == subpocket
                            and neighbor == frag2.ports[i].subpocket
                            and frag2.ports[i].bond_type == bond_type
                            and environment_match
                        ):
                            match.append([fragment_id1, fragment_id2])
                matches.append(
                    match
                )  # add valid matching pair to list of matching pairs
    return matches


def get_valid_pairs(fragment_library):

    """
    *copied and adapted from kinase_focused_fragment_library*
    Function preparing the fragment library to build pairs.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information

    Returns
    -------
    dict
        fragment library prepared for building valid pairs

    """
    data = {}  # (Fragments)
    frag_set = set()
    # only used in initialization for avoiding duplicates in fragment data set
    # (smiles & dummy atoms)

    # iterate through subpockets and fragments in subpockets
    # save subpocket_fragmentindex and dummy atoms, bonds etc
    for subpocket in fragment_library.keys():
        fragments = []
        for i, row in fragment_library[subpocket].iterrows():
            # get fragment and connecting subpockets
            fragment = row["ROMol_original"]
            fragment = Chem.RemoveHs(fragment)
            frag_id = f"{subpocket}_{i}"

            # store unique atom identifiers
            for a, atom in enumerate(fragment.GetAtoms()):
                frag_atom_id = f"{subpocket}_{a}"
                atom.SetProp("frag_atom_id", frag_atom_id)

            # get all dummy atoms of this fragment except the ones corresponding to the X pool
            dummy_atoms = [
                a
                for a in fragment.GetAtoms()
                if a.GetSymbol() == "*" and not a.GetProp("subpocket").startswith("X")
            ]
            if not dummy_atoms:
                continue

            frag_smiles, dummy_set = get_tuple(fragment, dummy_atoms)
            # check if this exact fragment has already been found
            if (frag_smiles, dummy_set) in frag_set:
                continue
            # if not, add this fragment to set of fragments
            frag_set.add((frag_smiles, dummy_set))

            # create dummy atom objects
            ports = [
                Port(
                    atom_id=dummy.GetProp("frag_atom_id"),
                    subpocket=subpocket,
                    neighboring_subpocket=dummy.GetProp("subpocket"),
                    bond_type=fragment.GetBondBetweenAtoms(
                        dummy.GetIdx(), dummy.GetNeighbors()[0].GetIdx()
                    ).GetBondType(),
                    environment=dummy.GetNeighbors()[0].GetProp("environment"),
                )
                for dummy in dummy_atoms
            ]

            # store fragment in constant data set
            fragment = Fragment(frag_id=frag_id, subpocket=subpocket, ports=ports)
            fragments.append(fragment)
        data[subpocket] = fragments

    n_frags = len(frag_set)

    print("Number of fragments: ", n_frags)

    return data


def get_tuple(fragment, dummy_atoms):

    """
    **copied from https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501
    For a given fragment, returns:
    - smiles string with generic dummy atoms (dummy labels removed)
    - dummy atoms as tuples of frag_atom_id and subpocket (of the dummy = neighboring subpocket of
    the fragment)
    Parameters
    ----------
    fragment: RDKit Mol object
    dummy_atoms: list(RDKit Atom objects)
        list of all dummy atoms of the fragment
    Returns
    -------
    String
        SMILES string of the fragment
    frozenset(tuple)
        frozenset of tuples for each dummy atom containing the frag_atom_id and the subpocket of
        the dummy
    """

    frag_smiles = fragment
    # replace dummys with generic dummys (without atom number)
    # dummy tuple: (frag_atom_id, neighboring_subpocket), e.g. (AP_4, FP)
    dummy_set = []
    for dummy in dummy_atoms:
        frag_smiles = Chem.ReplaceSubstructs(
            frag_smiles, Chem.MolFromSmiles(dummy.GetSmarts()), Chem.MolFromSmiles("*")
        )[0]
        dummy_tuple = dummy.GetProp("frag_atom_id"), dummy.GetProp("subpocket")
        dummy_set.append(dummy_tuple)
    frag_smiles = Chem.MolToSmiles(frag_smiles)

    dummy_set = frozenset(dummy_set)

    return frag_smiles, dummy_set


class Compound:

    """
    **copied from https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501
    Represents a combination of fragments including its dummy atoms
    Attributes
    ----------
    frag_ids: list(str)
        Strings representing the fragments that the molecule consists of
    subpockets: list(str)
        Subpockets that the molecule is targeting
    ports: list(Port)
        Port objects representing the dummy atoms of the molecule
    bonds: list(tuple(str))
        Bonds through which the fragments are connected.
        The bonds are stored as tuples of atom IDs.
    """

    def __init__(self, frag_ids, subpockets, ports, bonds):

        self.frag_ids = frag_ids
        self.subpockets = subpockets
        self.ports = ports
        self.bonds = bonds


class Fragment:

    """
    **copied from https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501
    Represents a single fragment from the fragment library
    Attributes
    ----------
    frag_id: str
        ID of the fragment: subpocket_ID, e.g. AP_5
    subpocket: str
        Subpocket that the fragment is targeting
    ports: list(Port)
        Port objects representing the dummy atoms of the fragment
    """

    def __init__(self, frag_id, subpocket, ports):

        self.frag_id = frag_id
        self.subpocket = subpocket  # list of targeted subpockets
        self.ports = ports  # list of Port objects


class Port:

    """
    **copied from https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501
    Represents a single dummy atom
    Attributes
    ----------
    atom_id: str
        frag_atom_id of the dummy atom
    subpocket: str
        Subpocket of the atom adjacent to the dummy atom (subpocket of the fragment containing
        the dummy)
    neighboring_subpocket: str
        Subpocket of the dummy atom
    bond_type: str
        Type of the bond connecting the dummy to its adjacent atom
    environment: str
        Type of the environment of the current fragment (of the adjacent atom)
    """

    def __init__(
        self, atom_id, subpocket, neighboring_subpocket, bond_type, environment
    ):

        self.atom_id = atom_id
        self.subpocket = subpocket
        self.neighboring_subpocket = neighboring_subpocket
        self.bond_type = bond_type
        self.environment = environment


class Combination:

    """
    **copied from https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501
    Comparable representation of a combination of fragments
    Attributes
    ----------
    frag_ids: frozenset(str)
        Strings representing the fragments that the molecule consists of
    bonds: frozenset(tuple(str))
        Bonds through which the fragments are connected.
        The bonds are stored as tuples of atom IDs.
    Methods
    ----------
    __eq__()
        Two Combination objects are equal if they consist of the same fragments which are
        connected through the same bonds.
    """

    def __init__(self, frag_ids, bonds=None):
        self.frag_ids = frag_ids
        self.bonds = bonds

    def __eq__(self, other):
        return self.frag_ids == other.frag_ids and self.bonds == other.bonds

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.frag_ids, self.bonds))
