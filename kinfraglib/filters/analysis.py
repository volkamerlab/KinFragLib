"""
Contains functions to analyze the results from filter steps
"""
import pandas as pd
from . import prefilters


def count_accepted_rejected(fragment_library, bool_column_name, filtername):
    """
    Function to count number of accepted and rejected fragments by one boolean column.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    bool_column_name : str
        string defining the column name where the bool values used for counting the accepted and
        rejected fragments are stored
    filtername : str
        string defining from which filter these bool column comes, e.g. "QED"

    Returns
    -------
    pandas.DataFrame
        number of accepted and number of rected fragments per subpocket
    """
    # seperate df into 0 and 1 in this column
    # concatenate fragment library and group fragments by the bool column
    fraglib_df = (
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby(bool_column_name, sort=False)
    )
    # count the df's grouped by subpocket
    accepted = pd.Series(
        fraglib_df.get_group(1).groupby("subpocket", sort=False).size(),
        name="accepted_" + filtername,
    )
    rejected = pd.Series(
        fraglib_df.get_group(0).groupby("subpocket", sort=False).size(),
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
    Function to count the number of accepted and rejected fragments.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    name : str
        string defining the column name what is counted

    Returns
    -------
    pandas.Series
        number of fragments per subpocket of the given dict
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
    Function to count number of fragments that are accepted by at least min_accepted filters
    defined in columns.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    min_accepted : int
        minimum number of accepted filters
    bool_column_name : str
        string defining the column names where the bool values are stored

    Returns
    -------
    dict
        containing a pandas.DataFrame including bool column if fragment is accepted by at least
        min_accepted filters which are defined in columns for each subpocket
    """
    fraglib_df = pd.concat(fragment_library).reset_index(drop=True)
    fraglib_df[name] = (fraglib_df.loc[:, columns].sum(axis=1) >= min_accepted).astype(int)
    return prefilters._make_df_dict(pd.DataFrame(fraglib_df))


def accepted_num_filters(fragment_library, colnames, filtername, max_num_accepted=1):
    """
    Function to count how many fragments are accepted by max_num_accepted or less filters.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets inculding all information
    colnames : list
        list containing strings with filter boolean column names
    filtername : str
        summarized filter name/ name of the resulting dataframe
    max_num_accepted : int
        maximum of accepted filters. By default max_num_accepted = 1

    Returns
    -------
    DataFrame
        Counting number of fragments that are accepted by max_num_accepted filters and less.
    """
    count_df = []
    cols_df = pd.concat(fragment_library).reset_index(drop=True)
    cols_df['sums'] = cols_df.loc[:, colnames].sum(axis=1)
    count_df.append(pd.Series(
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby("subpocket", sort=False)
        .size(),
        name="pre-filtered",
    ))
    for i in range(max_num_accepted, -1, -1):
        i_accepted = cols_df.loc[cols_df['sums'] == i]
        count_df.append(pd.Series(
            pd.concat(prefilters._make_df_dict(pd.DataFrame(i_accepted)))
            .reset_index(drop=True)
            .groupby("subpocket", sort=False)
            .size(),
            name=str("accepted by " + str(i)),
        ))
    counted_df = pd.concat(count_df, axis=1)
    counted_df = counted_df.fillna(0)
    counted_df = counted_df.astype(int)
    counted_df = counted_df.append(counted_df.sum().rename("Total"))
    counted_df = counted_df.style.set_caption(filtername)
    return(counted_df)
