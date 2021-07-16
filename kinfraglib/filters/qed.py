"""
Contains function to calculate the qed for each fragment.
"""
from rdkit import Chem
from . import check


def get_qed(fragment_library, cutoff_val=0.42, cutoff_crit=">"):
    """
    Calculates the Quantitavtive Estimate of Druglikeness.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    cutoff_val : str
        A value defining the cutoff for accepted/rejected fragments.
    cutoff_crit : str
        Defining wheter the QED value should be greater, equal, greater-equal, unequal, less or
        less equal compared to cutoff-val.

    Returns
    -------
    dict
        Containing
            A pandas.DataFrame with accepted fragments and their information.
            A pandas.DataFrame with rejected fragments and their information
            A dict containing a pandas.DataFrame for each subpocket with all fragments and an
            additional columns defining wether the fragment is accepted (1) or rejected (0)
            A list with the QED values
    """
    qedscores = []
    for pocket in fragment_library.keys():
        curqed = []
        for frag in fragment_library[pocket]["smiles"]:
            m = Chem.MolFromSmiles(frag)
            qed = Chem.QED.qed(m)
            curqed.append(qed)
        qedscores.append(curqed)
    qed_accepted, qed_rejected, fragment_library_bool = check.accepted_rejected(
        fragment_library,
        qedscores,
        cutoff_value=cutoff_val,
        cutoff_criteria=cutoff_crit,
        column_name="qed",
    )
    d = dict()
    d["qed_accepted"] = qed_accepted
    d["qed_rejected"] = qed_rejected
    d["fragment_library"] = fragment_library_bool
    d["qed"] = qedscores
    return d
