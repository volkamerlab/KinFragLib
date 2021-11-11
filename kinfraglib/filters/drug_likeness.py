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
    fragment_library : dict
        fragments organized in subpockets including all information
    min_fulfilled : int
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
    ro3_results = []        # array to fulfill the Ro3 results for each molecule
    num_fullfilled = []     # parameter to save the number of fulfilled Ro3 parameters
    for subpocket in fragment_library.keys():
        ro3_subpocket = []      # save ro3_booleans per subpocket
        for mol in fragment_library[subpocket]["ROMol"]:
            # call function to get for each Ro3 parameter a boolean if it is fulfilled
            ro3_subpocket.append(kfl_utils.get_ro3_from_mol(mol))
        ro3_results.append(ro3_subpocket)   # append Ro3 booleand to complete list
    # go through boolean list and save how many Ro3 parameters are filfilled
    for i in range(0, len(ro3_results)):
        num_sp = []
        for vals in ro3_results[i]:
            num_sp.append(sum(vals))
        num_fullfilled.append(num_sp)
    # add a boolean column to save whether the fragment gets accepted
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
    Calculates the Quantitative Estimate of Druglikeness.

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information
    cutoff_val : str
        A value defining the cutoff for accepted/rejected fragments. By default cutoff_val=0.42.
    cutoff_crit : str
        Defining whether the QED value should be ">", "<", ">=", "<=", "==" or "!=" compared to
        the cutoff-value. By default cutoff_crit=">".

    Returns
    -------
    dict
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional columns (bool_qed) defining whether the fragment is accepted (1) or rejected (0)
        and the calculated QED value for each fragment (qed).
    """
    qedscores = []      # save QED values
    # iterate through subpockets
    for subpocket in fragment_library.keys():
        curqed = []     # save current QED values from fragments in current subpocket
        for fragmol in fragment_library[subpocket]["ROMol"]:
            qed = Chem.QED.qed(fragmol)     # compute QED for each fragment
            curqed.append(qed)
        qedscores.append(curqed)       # save QED values
    # check if fragments accepted/rejected and add boolean column to fragment library
    fragment_library_bool = check.accepted_rejected(
        fragment_library,
        qedscores,
        cutoff_value=cutoff_val,
        cutoff_criteria=cutoff_crit,
        column_name="bool_qed",
    )
    # add column with QED values to the fragment library
    fragment_library_bool = utils.add_values(fragment_library_bool, qedscores, "qed")

    return fragment_library_bool
