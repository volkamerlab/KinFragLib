"""
Contains function to check which fragments are accepted or rejectes
"""

import pandas as pd
from . import prefilters
from . import building_blocks
import operator


def accepted_rejected(
    fragment_library,
    value_list,
    cutoff_value=0,
    cutoff_criteria="<",
    column_name="bool",
):
    """
    Go through values list and return a pandas.DataFrame of accepted/rejected fragments
    and a boolean list if fragment with this cutoff is rejected or accepted.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    value_list : list
        list of values calculated for filtering
    cutoff_value : int or float
        value defining the cutoff for accepting or rejecting a fragment
    cutoff_criteria : string of a basic operator
        defining if the rejected fragments values need to be >, <, >=, <=, == or != compared to
        the cutoff_value

    Returns
    -------
    pandas.DataFrames
        accepted/rejected ligands

    list of bools
        bool defining if this fragment is accepted or rejected
    """
    accepted = []
    rejected = []
    bools = []
    subpockets = list(fragment_library.keys())

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
        pocket = subpockets[i]
        # go through values in array
        for j in range(0, len(value_list[i])):
            val = value_list[i][j]
            # compare value with cutoff
            if ops[cutoff_criteria](val, cutoff_value):
                accepted.append(fragment_library[pocket].loc[j])
                bools.append(1)
            else:
                rejected.append(fragment_library[pocket].loc[j])
                bools.append(0)
    # save bool in fraglib df
    fragment_library_bool = building_blocks._add_bool_column(fragment_library, bools, column_name)
    return (
        prefilters._make_df_dict(pd.DataFrame(accepted)),
        prefilters._make_df_dict(pd.DataFrame(rejected)),
        fragment_library_bool,
    )
