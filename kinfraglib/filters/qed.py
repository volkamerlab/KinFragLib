"""
Contains function to calculate the qed for each fragment.
"""
from rdkit import Chem
from . import check
from . import utils


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
        Containing a pandas.DataFrame for each subpocket with all fragments and an
        additional columns defining wether the fragment is accepted (1) or rejected (0) and the
        calculated QED for each fragment.
    """
    qedscores = []
    for pocket in fragment_library.keys():
        curqed = []
        for frag in fragment_library[pocket]["smiles"]:
            m = Chem.MolFromSmiles(frag)
            qed = Chem.QED.qed(m)
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
