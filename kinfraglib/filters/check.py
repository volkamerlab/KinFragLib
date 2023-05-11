"""
Contains functions to check which fragments are accepted or rejected
"""
from . import synthesizability
import operator


def accepted_rejected(
    fragment_library,
    value_list,
    cutoff_value=0,
    cutoff_criteria="<",
    column_name="bool",
):
    """
    Go through value_list, compare it with the given cutoff and add a boolean column if
    fragments are accepted or rejected.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    value_list : list
        list of values calculated by a filtering step for filtering
    cutoff_value : int or float
        value defining the cutoff for accepting or rejecting a fragment
    cutoff_criteria : string of a basic operator
        defining if the rejected fragment values need to be ">", "<", ">=", "<=", "==" or "!="
        compared to the cutoff_value

    Returns
    -------
    dict
        fragment library containing a boolean column if the fragment was accepted (1) or
        rejected (0) by the filter according to the cutoff value.
    """
    bools = []

    # define operator list for comparison
    ops = {
        "<": operator.lt,
        ">": operator.gt,
        "==": operator.eq,
        "<=": operator.le,
        ">=": operator.ge,
        "!=": operator.ne,
    }
    # go through series indexes
    for i in range(0, len(value_list)):
        # go through values in array
        for j in range(0, len(value_list[i])):
            val = value_list[i][j]
            # compare value with cutoff
            if ops[cutoff_criteria](val, cutoff_value):
                bools.append(1)
            else:
                bools.append(0)
    # save bool column in in fragment library  df
    fragment_library_bool = synthesizability._add_bool_column(
        fragment_library, bools, column_name
    )
    return fragment_library_bool
