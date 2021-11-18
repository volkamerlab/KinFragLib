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
        fragments organized in subpockets including all information
    bool_column_name : str
        string defining the column name where the boolean values used for counting the accepted and
        rejected fragments are stored
    filtername : str
        string defining from which filter the given boolean column originates, e.g. "QED"

    Returns
    -------
    pandas.DataFrame
        number of accepted and number of rected fragments per subpocket
    """
    # concatenate fragment library and group fragments by the bool column
    fraglib_df = (
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby(bool_column_name, sort=False)
    )
    # count the fragments grouped by subpocket
    # accepted fragments
    accepted = pd.Series(
        fraglib_df.get_group(1).groupby("subpocket", sort=False).size(),
        name="accepted_" + filtername,
    )
    # rejected fragments
    rejected = pd.Series(
        fraglib_df.get_group(0).groupby("subpocket", sort=False).size(),
        name="rejected_" + filtername,
    )
    # remove NaN values and fill it with zeros
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

    # return the dataframe with the number of accepted and rejected fragments per subpockezt
    return accepted_rejected_df


def count_fragments(fragment_library, name="n_frags"):
    """
    Function to count the number of accepted and rejected fragments.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    name : str
        string defining the column name what is counted

    Returns
    -------
    pandas.Series
        number of fragments per subpocket of the given dict
    """
    # count and return the number of fragments per subpocket
    return pd.Series(
        pd.concat(fragment_library)
        .reset_index(drop=True)
        .groupby("subpocket", sort=False)
        .size(),
        name=name,
    )


def number_of_accepted(fragment_library, columns, min_accepted=1, name="bool"):
    """
    Function to count the number of fragments that are accepted by at least "min_accepted" filters
    specified in the given columns.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    columns : list
        of strings defining the column names where the boolean columns are stored
    min_accepted : int
        minimum number of accepted filters
    bool_column_name : str or list of str
        string or list of strings defining the column names where the boolean values are stored

    Returns
    -------
    dict
        containing a pandas.DataFrame including column where the number of accepted filters,
        defined in columns, is stored.
    """
    fraglib_df = pd.concat(fragment_library).reset_index(drop=True)
    # sum up number of accepted filters
    fraglib_df[name] = (fraglib_df.loc[:, columns].sum(axis=1) >= min_accepted).astype(
        int
    )
    return prefilters._make_df_dict(pd.DataFrame(fraglib_df))


def accepted_num_filters(fragment_library, colnames, filtername, max_num_accepted=1):
    """
    Function to count how many fragments are accepted by max_num_accepted or less filters.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    colnames : list
        list containing strings with the filter boolean column names
    filtername : str
        summarized filter name/ name of the resulting DataFrame
    max_num_accepted : int
        maximum of accepted filters. By default max_num_accepted = 1

    Returns
    -------
    DataFrame
        Counting number of fragments that are accepted by max_num_accepted filters and less.
    """
    # variable to store the number of accepted by num_accepted
    count_df = []
    # count by how many filters the fragments are accepted
    cols_df = pd.concat(fragment_library).reset_index(drop=True)
    cols_df["sums"] = cols_df.loc[:, colnames].sum(axis=1)
    # count the number of fragments in the pre-filtered library
    count_df.append(
        pd.Series(
            pd.concat(fragment_library)
            .reset_index(drop=True)
            .groupby("subpocket", sort=False)
            .size(),
            name="pre-filtered",
        )
    )
    # count number of accepted from max_num_accepted to zero
    for i in range(max_num_accepted, -1, -1):
        # get all fragments that are accepted by i filters
        i_accepted = cols_df.loc[cols_df["sums"] == i]
        # count them and store the numbers as a column in count_df
        count_df.append(
            pd.Series(
                pd.concat(prefilters._make_df_dict(pd.DataFrame(i_accepted)))
                .reset_index(drop=True)
                .groupby("subpocket", sort=False)
                .size(),
                name=str("accepted by " + str(i)),
            )
        )
    # create the final counted DataFrame
    counted_df = pd.concat(count_df, axis=1)
    # remove NA values and fill them with zeros
    counted_df = counted_df.fillna(0)
    # make numbers int not float (nicer to read)
    counted_df = counted_df.astype(int)
    # add a total number row at the end
    counted_df = counted_df.append(counted_df.sum().rename("Total"))
    # add the chosen title to the DataFrame
    counted_df = counted_df.style.set_caption(filtername)
    return counted_df
