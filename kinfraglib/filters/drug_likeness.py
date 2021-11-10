"""
Contains function to check for drug/lead likeness
"""
from rdkit import Chem
from kinfraglib import utils as kfl_utils
from . import check
from . import utils


def get_ro3_frags(fragment_library, min_fulfilled=6, cutoff_crit=">="):
    """
    Check if the fragments fulfill the rule of three criteria
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
        defined number of Ro3 parameters (bool_ro3).
    """
    ro3_results = []
    num_fullfilled = []
    all_fullfilled = []
    for subpocket in fragment_library.keys():
        ro3_subpocket = []
        for mol in fragment_library[subpocket]["ROMol"]:
            ro3_subpocket.append(kfl_utils.get_ro3_from_mol(mol))
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


def get_qed(fragment_library, cutoff_val=0.42, cutoff_crit=">"):
    """
    Calculates the Quantitavtive Estimate of Druglikeness.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets inculding all information
    cutoff_val : str
        A value defining the cutoff for accepted/rejected fragments.
    cutoff_crit : str
        Defining wheter the QED value should be greater, equal, greater-equal, unequal, less or
        less equal compared to cutoff-val.

    Returns
    -------
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional columns (bool_qed) defining whether the fragment is accepted (1) or rejected (0)
        and the calculated QED value for each fragment.
    """
    qedscores = []
    for pocket in fragment_library.keys():
        curqed = []
        for fragmol in fragment_library[pocket]["ROMol"]:
            qed = Chem.QED.qed(fragmol)
            curqed.append(qed)
        qedscores.append(curqed)
    fragment_library_bool = check.accepted_rejected(
        fragment_library,
        qedscores,
        cutoff_value=cutoff_val,
        cutoff_criteria=cutoff_crit,
        column_name="bool_qed",
    )
    fragment_library_bool = utils.add_values(fragment_library_bool, qedscores, "qed")

    return fragment_library_bool
