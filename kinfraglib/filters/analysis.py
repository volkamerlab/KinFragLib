"""
Contains functions to analyze the results from filter steps
"""
import pandas as pd
from . import prefilters
from . import pipeline_analysis
from kinfraglib import utils as kfl_utils
from IPython.display import display


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
        number of accepted and number of rejected fragments per subpocket
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
    Function to count how many fragments are accepted by max_num_accepted or fewer filters.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    colnames : list
        list containing strings with the filter boolean column names
    filtername : str
        summarized filter name/ name of the resulting DataFrame
    max_num_accepted : int
        maximum of accepted filters. By default, max_num_accepted = 1

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


def frag_in_subset(fragment_library_original, fragment_library_subset, colname):
    """
    Adding a boolean column to the fragment library if the fragments are contained in the fragment
    library subset.

    Parameters
    ----------
    fragment_libray : dict
        fragments organized in subpockets including all information
    fragment_libray : dict
        subset of the fragments organized in subpockets including all information
    colname : str
        name of the boolean column that will be created

    Returns
    -------
    dict
        fragments organized in subpocket with boolean column if the single fragments are included
        in the subset
    """
    # create dataframes from the fragment library dictionaries
    fragment_library_concat = pd.concat(fragment_library_original).reset_index(drop=True)
    fragment_library_reduced_concat = pd.concat(fragment_library_subset).reset_index(drop=True)

    bool_reduced = []       # variable to store the boolean columns
    # iterate through the fragment library
    for i, row in fragment_library_concat.iterrows():
        notfound = True
        # iterate through the fragment library subset
        for j, reduced_row in fragment_library_reduced_concat.iterrows():
            # compare the smiles, if they are equal fragment is in subset
            if row['smiles'] == reduced_row['smiles']:
                bool_reduced.append(1)
                notfound = False
                break
        if notfound:
            bool_reduced.append(0)
    # add the boolean column to the dataframe
    fragment_library_concat[colname] = bool_reduced
    # return the fragment library as dict
    fraglib = prefilters._make_df_dict(fragment_library_concat)
    return fraglib


def get_descriptors(fragment_library, fragment_library_reduced, fragment_library_custom):
    """
    Get #HBA #HBD, LogP and #Heavy Atoms for each fragment set and create a bar plot.

    Parameters
    ----------
    fragment_library : dict
        pre-filtered fragment library organized in subpockets
    fragment_library_reduced : dict
        reduced fragment library organized in subpockets
    fragment_library_custom : dict
        custom filtered fragment library organized in subpockets

    """
    descriptors = kfl_utils.get_descriptors_by_fragments(fragment_library)
    descriptors_median = descriptors.groupby('subpocket').median()
    descriptors_reduced = kfl_utils.get_descriptors_by_fragments(fragment_library_reduced)
    descriptors_reduced_median = descriptors_reduced.groupby('subpocket').median()
    descriptors_custom = kfl_utils.get_descriptors_by_fragments(fragment_library_custom)
    descriptors_custom_median = descriptors_custom.groupby('subpocket').median()

    all_descriptors = pd.concat(
        [
            descriptors_median,
            descriptors_reduced_median,
            descriptors_custom_median,
        ],
        axis=1,
        keys=["pre-filtered", "reduced", "custom"]
    )
    # style creates strange floats
    all_descriptors = all_descriptors.style.set_properties(
        **{"background-color": "lightgrey"},
        subset=["pre-filtered", "custom"],
    ).set_precision(precision=2)

    display(all_descriptors)

    print("\033[47;1m fragment library pre-filtered \033[0m")
    plt = pipeline_analysis.plot_fragment_descriptors(descriptors)
    # plt.title("fragment library pre-filtered")
    plt.show()

    print("\033[47;1m fragment  library reduced \033[0m")
    plt_reduced = pipeline_analysis.plot_fragment_descriptors(descriptors_reduced)
    # plt.title("fragment library reduced")
    plt_reduced.show()

    print("\033[47;1m fragment library custom \033[0m")
    plt_custom = pipeline_analysis.plot_fragment_descriptors(descriptors_custom)
    # plt.title("fragment library custom")
    plt_custom.show()


def get_descriptors_filters(fragment_library_filter_res, bool_keys):
    """
    Get #HBA #HBD, LogP and #Heavy Atoms for all fragments passing a filter step.

    Parameters
    ----------
    fragment_library_filter_res : dict
        pre-filtered fragment library organized in subpockets containing the filtering results
    bool_keys : list
        of strings containing the names of the boolean columns defining if a fragment passed a
        filter or not

    Returns
    ----------
    DataFrame
        containing the descriptors for each filtered set

    """
    # first calculate the descriptors from the pre-filtered library and plot them
    print("\033[47;1m pre-filtered \033[0m")
    descriptors = kfl_utils.get_descriptors_by_fragments(fragment_library_filter_res)
    descriptors_median = descriptors.groupby('subpocket').median()
    plt = pipeline_analysis.plot_fragment_descriptors(descriptors)
    plt.show()
    descriptor_dfs = {"pre-filtered": descriptors_median}   # add descriptors to a dataframe
    # iterate through the filters boolean columns, calculate the descriptor for passing fragments
    # and create the plots
    for bool_key in bool_keys:
        fraglib_concat = pd.concat(fragment_library_filter_res)
        fraglib_filter = fraglib_concat[fraglib_concat[bool_key] == 1]
        fraglib_filter = prefilters._make_df_dict(fraglib_filter)
        descriptors = kfl_utils.get_descriptors_by_fragments(fraglib_filter)
        descriptors_median = descriptors.groupby('subpocket').median()
        descriptor_dfs[bool_key] = descriptors_median   # add the descriptors to the descriptor df

        print("\033[47;1m " + bool_key + " \033[0m")
        plt = pipeline_analysis.plot_fragment_descriptors(descriptors)
        plt.show()
    # return the descriptors dataframe
    return(descriptor_dfs)


def filter_res_in_fraglib(fragment_library, filter_results):
    """
    Add the filtering results to the fragment library.

    Parameters
    ----------
    fragment_library : dict
        pre-filtered fragment library organized in subpockets containing the filtering results
    filter_results : DataFrame
        containing the fragments SMILES, the subpocket and the filtering results

    Returns
    ----------
    dict
        pre-filtered fragment library organized in subpockets containing the filtering results
    list
        of strings with the boolean column names for the filters

    """
    # set subpocket and smiles as index to add the filter results to the correct fragment
    filter_results = filter_results.set_index(["subpocket", "smiles"])

    fragment_library_concat = pd.concat(fragment_library)
    fragment_library_concat = fragment_library_concat.set_index(["subpocket", "smiles"])
    # merge the 2 dataframes
    fraglib_filters = fragment_library_concat.merge(
        filter_results,
        left_on=['subpocket', 'smiles'],
        right_on=['subpocket', 'smiles'],
        how="outer"
    )

    # get the list of boolean values, defining if a fragment is passing a specific filter or not
    frag_keys = fraglib_filters.keys()
    frag_keys.to_list()
    bool_keys = [x for x in frag_keys if "bool" in x]

    # set index to subpocket again and crate dict
    fraglib_filters = fraglib_filters.reset_index()
    fraglib_filters.set_index(["subpocket"])

    fragment_library_filter_res = prefilters._make_df_dict(fraglib_filters)

    # return fragment library dict with the filtering results and the boolean keys
    return fragment_library_filter_res, bool_keys
