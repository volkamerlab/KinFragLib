"""
Contains functions to analyse the results from filter steps
"""
import pandas as pd
from . import prefilters


def count_accepted_rejected(fragment_library, bool_column_name, filtername):
    """
    Function to count number of accepted and rejected fragments by boolean column.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    bool_column_name : str
        string defining the column name where the bool values are stored
    filtername : str
        string defining from which filter these bool column comes

    Returns
    -------
    pandas.DataFrame
        number of accepted and number of rected fragments per subpocket
    """
    # seperate df into 0 and 1 in this column
    df = (
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby(bool_column_name, sort=False)
    )
    # count the df's grouped by subpocket
    accepted = pd.Series(
        df.get_group(1).groupby("subpocket", sort=False).size(),
        name="accepted_" + filtername,
    )
    rejected = pd.Series(
        df.get_group(0).groupby("subpocket", sort=False).size(),
        name="rejected_" + filtername,
    )
    # remove NaN values
    accepted_rejected_df = pd.concat([accepted, rejected], axis=1)
    accepted_rejected_df[str("rejected_" + filtername)] = accepted_rejected_df[
        str("rejected_" + filtername)
    ].fillna(0)
    accepted_rejected_df[str("rejected_" + filtername)] = accepted_rejected_df[
        str("rejected_" + filtername)
    ].astype(int)
    accepted_rejected_df[str("accepted_" + filtername)] = accepted_rejected_df[
        str("accepted_" + filtername)
    ].fillna(0)
    accepted_rejected_df[str("accepted_" + filtername)] = accepted_rejected_df[
        str("accepted_" + filtername)
    ].astype(int)
    return accepted_rejected_df


def count_fragments(fragment_library, name="n_frags"):
    """
    Function to count number of accepted and rejected fragments by boolean column.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    name : str
        string defining the column name where the bool values are stored

    Returns
    -------
    pandas.Series
        number of accepted and number of rected fragments per subpocket
    """
    return pd.Series(
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby("subpocket", sort=False)
        .size(),
        name=name,
    )


def number_of_accepted(fragment_library, columns, min_accepted=1, name="bool"):
    """
    Function to count number of accepted and rejected fragments by boolean column.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    min_accepted : int
        minimum number of accepted filters
    bool_column_name : str
        string defining the column name where the bool values are stored

    Returns
    -------
    dict
        containing a pandas.DataFrame including bool column for each subpocket
    """
    df = pd.concat(fragment_library).reset_index(drop=True)
    df[name] = (df.loc[:, columns].sum(axis=1) >= min_accepted).astype(int)
    return prefilters._make_df_dict(pd.DataFrame(df))
