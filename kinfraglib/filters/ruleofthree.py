"""
Contains function to check the rule of three parameters
"""
from rdkit import Chem
from kinfraglib import utils
from . import check


def get_ro3_frags(fragment_library, min_fulfilled=6, cutoff_crit=">="):
    """
    Check the rule of three parameters
        - molecular weight <300
        - logp <=3
        - number of hydrogen bond acceptors <=3
        - number of hydrogen bond donors <=3
        - number of rotatable bonds <=3
        - polar surface area <= 60

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    min_fiulfilled : int
        defining the minimum number of Rule of Three Criteria that need to be fulfilled to be
        accepted. Default: min_fulfilled=6
    cutoff_crit : str
        Per default ">="

    Returns
    -------
    dict
        Containing
            A pandas.DataFrame with accepted fragments and their information.
            A pandas.DataFrame with rejected fragments and their information
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional columns defining wether the fragment is accepted (1) or rejected (0)
            A list with the Ro3 results
    """
    ro3_results = []
    num_fullfilled = []
    all_fullfilled = []
    for subpocket in fragment_library.keys():
        ro3_subpocket = []
        for smiles in fragment_library[subpocket]["smiles"]:
            m = Chem.MolFromSmiles(smiles)
            ro3_subpocket.append(utils.get_ro3_from_mol(m))
        ro3_results.append(ro3_subpocket)
    for i in range(0, len(ro3_results)):
        num_sp = []
        num_bools = []
        for vals in ro3_results[i]:
            num_sp.append(sum(vals))
            if sum(vals) == 6:
                num_bools.append(1)
            else:
                num_bools.append(0)
        all_fullfilled.append(num_bools)
        num_fullfilled.append(num_sp)

    ro3_accepted, ro3_rejected, fragment_library_bool = check.accepted_rejected(
        fragment_library,
        num_fullfilled,
        cutoff_value=min_fulfilled,
        cutoff_criteria=cutoff_crit,
        column_name="ro3",
    )
    d = dict()
    d["ro3_accepted"] = ro3_accepted
    d["ro3_rejected"] = ro3_rejected
    d["fragment_library"] = fragment_library_bool
    d["ro3"] = ro3_results
    return d
