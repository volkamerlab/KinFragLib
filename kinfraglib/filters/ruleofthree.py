"""
Contains function to check the rule of three parameters
"""
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
        accepted. By default min_fulfilled=6.
    cutoff_crit : str
        Cutoff criterium, defining if the number of fulfilled parameters is ">", "<", "==", ">="
        or "<=" than min_fulfilled. By default cutoff_crit=">=".

    Returns
    -------
    dict
        fragment library organized in subpockets containing a boolean column if they fulfill the
        defined number of Ro3 parameters.
    """
    ro3_results = []
    num_fullfilled = []
    all_fullfilled = []
    for subpocket in fragment_library.keys():
        ro3_subpocket = []
        for mol in fragment_library[subpocket]["ROMol"]:
            ro3_subpocket.append(utils.get_ro3_from_mol(mol))
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

    fragment_library_bool = check.accepted_rejected(
        fragment_library,
        num_fullfilled,
        cutoff_value=min_fulfilled,
        cutoff_criteria=cutoff_crit,
        column_name="bool_ro3",
    )
    return fragment_library_bool
