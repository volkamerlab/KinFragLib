"""
Contains plotting functionalities for different filters.
"""

import statistics
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, MACCSkeys
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from IPython.display import display
from kinfraglib import utils as kfl_utils
from collections import Counter


def make_hists(
    fragment_library, colname, filtername=None, plot_stats=True, cutoff=None
):
    """
    Creates a histogram for each subpocket for the given values.

    Parameters
    ----------
    fragment_library : dict
        fragment library organized in subpockets.

    colname : str
        Name of the column where values for creating the histograms are stored.

    filtername : str
        name of the filter used as title creating the values plotted. By default, filtername=None,
        meaning no title will be displayed.

    plot_stats : boolean
        defining if a box with min, max and mean value will be displayed in each plot. By default
        plot_stats=True.

    cutoff : int or float
        cutoff value for drawing a cutoff line to the plots. By default, cutoff=None, meaning no
        cutoff line is plotted.
    """
    # get even number if number of plots not even
    num_plots = round(len(fragment_library.keys()) + 0.5)
    plt.figure(figsize=(20, 22))
    # create grid to place plots next to each other
    gs = gridspec.GridSpec(int(num_plots / 2), int(num_plots / 2))
    # save keys
    keys = list(fragment_library.keys())
    subpocket_num = 0
    # create one plot for each subpocket
    for i in range(0, 2):
        for j in range(0, int((num_plots) / 2)):
            if (i * 4) + j <= num_plots and subpocket_num < len(
                fragment_library.keys()
            ):
                ax = plt.subplot(gs[i, j])
                ax.hist(
                    fragment_library[keys[subpocket_num]][colname],
                    facecolor="#04D8B2",
                    edgecolor="#808080",
                )
                ax.set_title(keys[subpocket_num])
                if (
                    plot_stats
                ):  # add statistics box (max, min, mean value per subpocket)
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="mean: "
                        + str(  # noqa: W504
                            round(
                                statistics.mean(
                                    fragment_library[keys[subpocket_num]][colname]
                                ),
                                ndigits=2,
                            )
                        ),  # noqa: E501
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="min: "
                        + str(  # noqa: W504
                            round(min(fragment_library[keys[subpocket_num]][colname]))
                        ),
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="max: "
                        + str(  # noqa: W504
                            round(max(fragment_library[keys[subpocket_num]][colname]))
                        ),
                    )
                    plt.legend()
                if cutoff is not None:  # if a cutoff is given draw a red line
                    plt.axvline(x=cutoff, color="r", linestyle="-")
                if filtername is not None:
                    plt.xlabel(filtername)
                plt.ylabel("Number of fragments")  # set yaxis label
                subpocket_num = subpocket_num + 1  # go to next subpocket
    plt.suptitle(filtername)  # set filtername as title over the plots
    plt.show()  # show the plots, needed when called in function other not in notebook


def make_retro_hists(
    fragment_library, colname, filtername=None, plot_stats=True, cutoff=None
):
    """
    Creates a histogram for each subpocket for defined values.
    Parameters
    ----------
    fragment_library : dict
        fragment library organized in subpockets
    colname : str
        Name of the column where values for creating histograms are stored
    filtername : str
        name of the filter used as title creating the values plotted
    cutoff : int or float
        cutoff value for drawing a cutoff line to the plots
    """
    # only consider subpockets that are not empty after filtering
    fragment_library = {
        sp: data
        for sp, data in fragment_library.items()
        if any(fragment_library[sp][colname])
    }

    # get even number if number of plots not even
    num_plots = round(len(fragment_library.keys()) + 0.5)
    plt.figure(figsize=(25, 29))
    gs = gridspec.GridSpec(int(num_plots / 2), int(num_plots / 2))
    keys = list(fragment_library.keys())
    subpocket_num = 0
    for i in range(0, 2):
        for j in range(0, int((num_plots) / 2)):
            if (i * 4) + j <= num_plots and subpocket_num < len(
                fragment_library.keys()
            ):
                cur_data = fragment_library[keys[subpocket_num]][colname]
                cur_binsize = round(max(cur_data) / 9)
                bin_lst = list(range(0, max(cur_data) + cur_binsize, cur_binsize))
                bin_lst.pop(0)
                bin_lst = [-(cur_binsize), 0.1] + bin_lst
                medians = []
                bin_label = []
                for x in range(0, len(bin_lst) - 1):
                    if x == 0:
                        bin_str = "[0]"
                    elif x == 1:
                        bin_str = "(%s, %s)" % (0, bin_lst[x + 1])
                    elif x == len(bin_lst) - 1:
                        bin_str = "[%s, %s]" % (bin_lst[x], bin_lst[x + 1])
                    else:
                        bin_str = "[%s, %s)" % (bin_lst[x], bin_lst[x + 1])
                    bin_label.append(bin_str)
                    cur_bins = [bin_lst[x], bin_lst[x + 1]]
                    medians.append(statistics.median(cur_bins))
                medians = [round(num) for num in medians]
                label_pos = []
                for label in medians:
                    label_pos.append(label)
                ax = plt.subplot(gs[i, j])
                N, _, patches = ax.hist(cur_data, bins=bin_lst, rwidth=0.9)
                ax.set_xticks(label_pos)
                ax.set_xticklabels(bin_label, rotation=90)
                ax.set_title(keys[subpocket_num])
                patches[0].set_facecolor("r")
                for count, patch in zip(N, patches):
                    ax.annotate(
                        str(int(count)),
                        xy=(patch.get_x() + (cur_binsize / 2) - 1, patch.get_height()),
                        ha="center",
                        va="bottom",
                    )
                if plot_stats:
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="mean: "
                        + str(  # noqa: W504
                            round(
                                statistics.mean(
                                    fragment_library[keys[subpocket_num]][colname]
                                )
                            )
                        ),  # noqa: E501
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="min: "
                        + str(  # noqa: W504
                            round(min(fragment_library[keys[subpocket_num]][colname]))
                        ),
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="max: "
                        + str(  # noqa: W504
                            round(max(fragment_library[keys[subpocket_num]][colname]))
                        ),
                    )
                    plt.legend()
                if filtername is not None:
                    plt.xlabel(filtername)
                plt.ylabel("Number of fragments")
                plt.xlabel("Number of retrosynthetic routes")
                subpocket_num = subpocket_num + 1
    plt.suptitle(filtername)
    plt.show()


def retro_routes_fragments(fragment_library, evaluate, subpocket, molsPerRow=10):
    """
    Creates an image of the fragments for the given subpocket
    a) without a retrosynthetic route found if evaluate="none"
    b) max. 10 fragments with the most retrosynthetic routes found.
    ----------
    fragment_library : dict
        fragment library organized in subpockets

    evaluate : str
        "none" or "max", defining if the fragments without retrosynthetic route ("none") or
        the fragments with the most retrosynthetic routes found ("max") will be shown

    subpocket : str
        defining the fragments from which subpocket will be shown

    molsPerRow : int
        defining how many molecules are displayed in one row. By default, molsPerRow=10
    """
    if evaluate == "none":
        num_fragments = len(
            pd.Series(
                fragment_library[subpocket][
                    fragment_library[subpocket]["retro_count"] == 0
                ].ROMol
            )
        )
        print(
            "%s %s fragments with no retrosynthetic route found"
            % (num_fragments, subpocket)
        )
        if num_fragments == 0:
            return
        img = Draw.MolsToGridImage(
            pd.Series(
                fragment_library[subpocket][
                    fragment_library[subpocket]["retro_count"] == 0
                ].ROMol
            ),
            molsPerRow=molsPerRow,
            maxMols=len(
                pd.Series(
                    fragment_library[subpocket][
                        fragment_library[subpocket]["retro_count"] == 0
                    ].ROMol
                ),
            ),
            returnPNG=False,
        )
        return img
    elif evaluate == "max":
        num_fragments = len(
            pd.Series(
                fragment_library[subpocket][
                    fragment_library[subpocket]["retro_count"] > 0
                ][0:10].ROMol
            )
        )
        print(
            "%s %s fragments with the most retrosynthetic routes found"
            % (
                num_fragments,
                subpocket,
            )
        )
        if num_fragments == 0:
            return
        print("legend: number of retrosynthetic routes found")
        img = Draw.MolsToGridImage(
            pd.Series(
                fragment_library[subpocket][
                    fragment_library[subpocket]["retro_count"] > 0
                ]
                .sort_values("retro_count", ascending=False)[0:10]
                .ROMol
            ),
            molsPerRow=molsPerRow,
            legends=list(
                fragment_library[subpocket]
                .sort_values("retro_count", ascending=False, ignore_index=True)[0:10][
                    "retro_count"
                ]
                .astype(str)
            ),
            returnPNG=False,
        )
        return img


def sample_subset(fragment_library, sample):
    """
    Samples a random subset of the given fragment library
    ----------
    fragment_library : dict
        fragment library organized in subpockets
    sample : float
        fraction of data to be sampled

    """

    fragment_library_subset = {}
    if sample != 1.0:
        for subpocket in fragment_library.keys():
            sample_num = int(len(fragment_library[subpocket]) * sample)
            fragment_library_subset[subpocket] = fragment_library[subpocket].sample(
                sample_num, random_state=1
            )
    return fragment_library_subset


def create_tsne_embeddings(fragment_library):
    """
    Creates the t-SNE embedding for all following t-SNE plots
    ----------
    fragment_library : dict
        fragment library organized in subpockets

    """

    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    fragment_library_concat["maccs"] = fragment_library_concat.ROMol.apply(
        MACCSkeys.GenMACCSKeys
    )

    pca = PCA(n_components=30)
    crds = pca.fit_transform(list(fragment_library_concat["maccs"]))

    crds_embedded = TSNE(
        n_components=2, init="pca", learning_rate="auto"
    ).fit_transform(crds)
    return crds_embedded


def create_tsne_plots(crds_embedded, fragment_library):
    """
    Creates t-SNE plots comparing
    a) pre-filtered and reduced fragment library
    b) pre-filtered and custom filtered fragment library
    c) pre-filtered, reduced and custom fragment library

    and prints number of fragments in the subsets.
    ----------
    crds_embedded : list
        t-SNE embedding of the fragment library
    fragment_library : dict
        fragment library organized in subpockets containing boolean columuns `bool_reduced`and
        `bool_custom`defining if the fragments are part of the subsets

    """
    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    tsne_df = pd.DataFrame(crds_embedded, columns=["X", "Y"])
    # add bool column from filtering steps here
    tsne_df["reduced"] = fragment_library_concat["bool_reduced"]
    tsne_df["custom"] = fragment_library_concat["bool_custom"]
    # create column defining if fragment is
    # *excluded in both subsets (0)
    # *included in custom (1)
    # *included in reduced (2)
    # *included in both subsets (3)
    bool_compare = []
    for i, row in fragment_library_concat.iterrows():
        if row["bool_reduced"] == 0 and row["bool_custom"] == 0:
            bool_compare.append(0)
        elif row["bool_reduced"] == 0 and row["bool_custom"] == 1:
            bool_compare.append(1)
        elif row["bool_reduced"] == 1 and row["bool_custom"] == 0:
            bool_compare.append(2)
        elif row["bool_reduced"] == 1 and row["bool_custom"] == 1:
            bool_compare.append(3)
    tsne_df["compare"] = bool_compare
    num0 = len(tsne_df[tsne_df["compare"] == 0])
    num1 = len(tsne_df[tsne_df["compare"] == 1])
    num2 = len(tsne_df[tsne_df["compare"] == 2])
    num3 = len(tsne_df[tsne_df["compare"] == 3])

    # create tsne plots
    plt.figure(figsize=(18, 10))
    plt.subplot(2, 2, 1)
    sns.scatterplot(
        data=tsne_df.query("reduced == 0"),
        x="X",
        y="Y",
        color="lightcoral",
        alpha=0.5,
        label="excluded",
    ).set_title("pre_filtered vs. reduced")
    sns.scatterplot(
        data=tsne_df.query("reduced == 1"),
        x="X",
        y="Y",
        color="green",
        alpha=0.5,
        label="included",
    )

    plt.subplot(2, 2, 2)
    sns.scatterplot(
        data=tsne_df.query("custom == 0"),
        x="X",
        y="Y",
        color="lightcoral",
        alpha=0.5,
        label="excluded",
    ).set_title("pre-filtered vs. custom")
    sns.scatterplot(
        data=tsne_df.query("custom == 1"),
        x="X",
        y="Y",
        color="green",
        alpha=0.5,
        label="included",
    )

    plt.subplot(2, 2, 3)
    sns.scatterplot(
        data=tsne_df.query("compare == 0"),
        x="X",
        y="Y",
        color="lightcoral",
        alpha=0.5,
        label="excluded in both subsets",
    ).set_title("pre-filtered vs. reduced vs. custom")
    sns.scatterplot(
        data=tsne_df.query("compare == 1"),
        x="X",
        y="Y",
        color="orange",
        alpha=0.5,
        label="included in custom subset",
    )
    sns.scatterplot(
        data=tsne_df.query("compare == 2"),
        x="X",
        y="Y",
        color="lightblue",
        alpha=0.5,
        label="included in reduced subset",
    )
    sns.scatterplot(
        data=tsne_df.query("compare == 3"),
        x="X",
        y="Y",
        color="green",
        alpha=0.5,
        label="included in both subsets",
    )
    plt.legend(loc="upper right", bbox_to_anchor=(1.425, 1), ncol=1)

    plt.show()
    num_lists = (len(tsne_df["compare"]), num0, num1, num2, num3)
    print(
        """%s Pre-filtered fragments.
        Number of fragments excluded in both datasets: %s
        Number of fragments included in the custom dataset and excluded in the reduced dataset: %s
        Number of fragments included in the reduced dataset and excluded in the custom dataset: %s
        Number of fragments in both datasets: %s """
        % (num_lists)
    )
    tsne_df["smiles"] = fragment_library_concat["smiles"]
    return tsne_df


def create_tsne_plots_filters(crds_embedded, fragment_library, saved_filter_results):
    """
    Creates t-SNE plots with accepted (green) and rejected (red) fragments for each filtering step.

    ----------
    crds_embedded : list
        t-SNE embeddings of the fragment library
    fragment_library : dict
        fragment library organized in subpockets containing boolean columuns
    saved_filter_results : dataframe
        loaded file with saved filter results

    """
    # get sampled subset from saved filter results
    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    bool_samples = [
        i
        for i, smi in enumerate(saved_filter_results["smiles"])
        if smi in list(fragment_library_concat["smiles"])
    ]
    saved_filter_results = saved_filter_results.iloc[bool_samples]

    tsne_df = pd.DataFrame(crds_embedded, columns=["X", "Y"])
    # add bool column from filter steps
    filters = []
    if "bool_pains" in saved_filter_results.columns:
        tsne_df["pains"] = saved_filter_results["bool_pains"]
        filters.append("pains")
    if "bool_brenk" in saved_filter_results.columns:
        tsne_df["brenk"] = saved_filter_results["bool_brenk"]
        filters.append("brenk")
    if "bool_ro3" in saved_filter_results.columns:
        tsne_df["ro3"] = saved_filter_results["bool_ro3"]
        filters.append("ro3")
    if "bool_qed" in saved_filter_results.columns:
        tsne_df["qed"] = saved_filter_results["bool_qed"]
        filters.append("qed")
    if "bool_bb" in saved_filter_results.columns:
        tsne_df["bb"] = saved_filter_results["bool_bb"]
        filters.append("bb")
    if "bool_syba" in saved_filter_results.columns:
        tsne_df["syba"] = saved_filter_results["bool_syba"]
        filters.append("syba")
    if "bool_retro" in saved_filter_results.columns:
        tsne_df["retro"] = saved_filter_results["bool_retro"]
        filters.append("retro")

    tsne_df["smiles"] = saved_filter_results["smiles"]

    # create the plots for all filters
    plt.figure(figsize=(15, 21))
    i = 0
    for filter in filters:
        i = i + 1
        plt.subplot(4, 2, i)
        sns.scatterplot(
            data=tsne_df.query("%s == 1" % filter),
            x="X",
            y="Y",
            color="green",
            alpha=0.5,
            label="accepted",
        ).set_title(filter)
        sns.scatterplot(
            data=tsne_df.query("%s == 0" % filter),
            x="X",
            y="Y",
            color="lightcoral",
            label="rejected",
        )
        # plt.show()
    return tsne_df


def connection_frequencies(
    fragment_library, fragment_library_reduced, fragment_library_custom
):
    """
    Calculates the connection frequencies between the subpockets for all three subsets and
    creates a barplot.

    *part of the code copied from https://github.com/sonjaleo/KinFragLib/blob/master/notebooks/2_2_fragment_analysis_statistics.ipynb   # noqa: E501

    ----------
    fragment_library : dict
        pre-filtered fragment library organized in subpockets
    fragment_library_reduced : dict
        reduced fragment library organized in subpockets
    fragment_library_custom : dict
        custom filtered fragment library organized in subpockets

    Returns
    ---------
    dataframe
        with the connection frequencies for every existing connection between subpockets for all
        three subsets

    """
    # fragment library pre-filtered
    fragment_library_concat = pd.concat(fragment_library)
    connections_by_fragment = kfl_utils.get_connections_by_fragment(
        fragment_library_concat
    )
    connections_by_ligand = connections_by_fragment.groupby(
        ["kinase", "complex_pdb", "ligand_pdb"]
    )["connections_name"].sum()
    connections_by_ligand_count = connections_by_ligand.apply(lambda x: Counter(x))
    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count = pd.Series(
        Counter(connections_by_ligand_count.apply(list).sum())
    )
    connections_across_ligands_count.name = "count_pre-filtered"

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency = connections_across_ligands_count.apply(
        lambda x: round((x / connections_by_ligand_count.shape[0] * 100), 1)
    )
    connections_across_ligands_frequency.name = "frequency_pre-filtered"

    # Concatenate count and frequency data to DataFrame
    connections_across_ligands = pd.concat(
        [connections_across_ligands_count, connections_across_ligands_frequency],
        axis=1,
    )

    # fragment library reduced
    fragment_library_reduced_concat = pd.concat(fragment_library_reduced)
    connections_by_fragment_reduced = kfl_utils.get_connections_by_fragment(
        fragment_library_reduced_concat
    )

    connections_by_ligand_reduced = connections_by_fragment_reduced.groupby(
        ["kinase", "complex_pdb", "ligand_pdb"]
    )["connections_name"].sum()
    connections_by_ligand_count_reduced = connections_by_ligand_reduced.apply(
        lambda x: Counter(x)
    )

    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count_reduced = pd.Series(
        Counter(connections_by_ligand_count_reduced.apply(list).sum())
    )
    connections_across_ligands_count_reduced.name = "count_reduced"

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency_reduced = (
        connections_across_ligands_count_reduced.apply(
            lambda x: round((x / connections_by_ligand_count_reduced.shape[0] * 100), 1)
        )
    )
    connections_across_ligands_frequency_reduced.name = "frequency_reduced"

    # Concatenate count and frequency data to DataFrame
    connections_across_ligands_reduced = pd.concat(
        [
            connections_across_ligands_count_reduced,
            connections_across_ligands_frequency_reduced,
        ],
        axis=1,
    )

    # fragment library custom filtered
    fragment_library_custom_concat = pd.concat(fragment_library_custom)
    connections_by_fragment_custom = kfl_utils.get_connections_by_fragment(
        fragment_library_custom_concat
    )

    connections_by_ligand_custom = connections_by_fragment_custom.groupby(
        ["kinase", "complex_pdb", "ligand_pdb"]
    )["connections_name"].sum()
    connections_by_ligand_count_custom = connections_by_ligand_custom.apply(
        lambda x: Counter(x)
    )

    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count_custom = pd.Series(
        Counter(connections_by_ligand_count_custom.apply(list).sum())
    )
    connections_across_ligands_count_custom.name = "count_custom-filtered"

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency_custom = (
        connections_across_ligands_count_custom.apply(
            lambda x: round((x / connections_by_ligand_count_custom.shape[0] * 100), 1)
        )
    )
    connections_across_ligands_frequency_custom.name = "frequency_custom_filtered"

    # Concatenate count and frequency data to DataFrame
    connections_across_ligands_custom = pd.concat(
        [
            connections_across_ligands_count_custom,
            connections_across_ligands_frequency_custom,
        ],
        axis=1,
    )

    # combine connection frequencies in one dataframe
    frequencies = pd.concat(
        [
            connections_across_ligands["frequency_pre-filtered"],
            connections_across_ligands_reduced["frequency_reduced"],
            connections_across_ligands_custom["frequency_custom_filtered"],
        ],
        axis=1,
    )

    # create the barplot to show the connection frequencies for each subset and each connection
    ax = frequencies.plot.bar()
    fig = ax.get_figure()

    fig.set_figheight(5)
    fig.set_figwidth(13)

    ax.set_xlabel("Subpocket")
    ax.set_ylabel("Connection Frequency")
    ax.set_title("Connection Frequencies of the different subsets")

    res = pd.concat(
        [
            connections_across_ligands,
            connections_across_ligands_reduced,
            connections_across_ligands_custom,
        ],
        axis=1,
    )

    res = res.fillna(0)
    res["count_custom-filtered"] = res["count_custom-filtered"].astype(int)
    return res


def num_frags_development(filter_res):
    """
    Count the number of fragments passing each custom filter step

    ----------
    filter_res : dataframe
        Contains the calculated values and the boolean for each filtering step if a fragment was
        accepted or not.

    Returns
    ---------
    dataframe
        with the number of fragments per subpocket after each filtering step

    """
    # get the column names
    frag_keys = filter_res.keys()
    frag_keys.to_list()
    # keep only the boolean column names (we do not need the computed values here)
    bool_keys = [x for x in frag_keys if "bool" in x]
    # create a dataframe to store the number of fragments left after each filtering step
    update_results = pd.DataFrame()
    # add number of fragments for the pre-filtered subset we are starting with
    update_results["pre-filtered"] = (
        filter_res.reset_index().groupby("subpocket", sort=False).size()
    )
    # go through all boolean columns after one another and count the number of fragments passing
    for bool_key in bool_keys:
        filter_res = filter_res.loc[filter_res[bool_key] == 1]
        update_results[bool_key] = (
            filter_res.reset_index().groupby("subpocket", sort=False).size()
        )

    # create a bar plot showing the numbers of fragments passing
    ax = update_results.plot.bar(width=0.9)
    fig = ax.get_figure()

    fig.set_figheight(7)
    fig.set_figwidth(13)

    ax.set_xlabel("Subpocket")
    ax.set_ylabel("Number of fragments")
    ax.set_title(
        "Development of the number of fragments per subpocket after each filter step"
    )
    # return dataframe with number of fragments after each filtering step.
    return update_results


SUBPOCKET_COLORS = {
    "AP": "purple",
    "FP": "forestgreen",
    "SE": "c",
    "GA": "tab:orange",
    "B1": "tab:blue",
    "B2": "darkslateblue",
    "X": "grey",
}


def draw_clusters(clustered_fragments, show_subpockets=False):
    # function copied from https://github.com/volkamerlab/KinFragLib/blob/master/notebooks/2_3_fragment_analysis_most_common_fragments.ipynb  # noqa: E501
    """
    Draw fragments sorted by descending cluster size and fragment count.

    Parameters
    ----------
    most_common_fragments : pandas.DataFrame
        Most common fragments (ID, SMILES, ROMol, cluster ID, fragment count).
    show_subpockets : boolean
        If show_subpockets=True the subpockets that contain this fragment are shown. By default,
        show_subpockets=False.

    Returns
    -------
    PIL.PngImagePlugin.PngImageFile
        Image of fragments sorted by descending cluster size.
    """
    # clustered_fragments = clustered_fragments.sort_values("subpocket_count", ascending=False)
    if show_subpockets:
        legend = [
            f"{row.cluster_id} | {row.subpocket_count} |  {row.subpockets}"
            for index, row in clustered_fragments.iterrows()
        ]
    else:
        legend = [
            f"{row.cluster_id} | {row.subpocket_count}"
            for index, row in clustered_fragments.iterrows()
        ]
    img = Draw.MolsToGridImage(
        list(clustered_fragments.ROMol),
        legends=legend,
        molsPerRow=7,
        maxMols=100,
        subImgSize=(170, 170),
        useSVG=True,
    )
    if show_subpockets:
        print("Legend: cluster ID | subpocket count | subpockets")
    else:
        print("Legend: cluster ID | subpocket count")

    return img


def plot_fragment_descriptors(descriptors, ylims):
    """
    Plot fragment descriptors.

    Parameters
    ----------
    descriptors : dataframe
        descriptors for each fragment

    Returns
    ---------
        Plot of the descriptors, colored by subpocket.
    """
    # copied from utils without saving img

    plt.figure(figsize=(25, 6))

    for i, descriptor_name in enumerate(descriptors.columns[3:]):

        plt.subplot(1, 4, i + 1)
        ax = sns.boxplot(
            x="subpocket",
            y=descriptor_name,
            hue="subpocket",
            legend=False,
            data=descriptors,
            palette=SUBPOCKET_COLORS,
            medianprops={"linewidth": 3, "linestyle": "-"},
        )
        ax.set_ylim((ylims[i][0], ylims[i][1]))
        plt.ylabel(descriptor_name, fontsize=16)
        plt.xlabel("Subpocket", fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

    return plt


def plot_fragment_similarity(similarities_by_groups, library_names, group_name):
    """
    Plot fragment similarity by category, such as subpocket or kinase group.

    Parameters
    ----------
    similiarities_by_groups : dataframe
        fragment similarity per subpocket/kinase group

    similiarities_by_groups : list
        containing strings with library names used as titles

    group_name : str
        the group that is plotted
    """
    # copied from kinfraglib/utils.py and modified

    plt.figure(figsize=(20, 5))
    num_plots = len(similarities_by_groups)
    i = 0
    for similarities_by_group in similarities_by_groups:
        plt.subplot(1, num_plots, i + 1)
        try:
            sns.boxplot(
                x=similarities_by_group.columns[1],
                y=similarities_by_group.columns[0],
                hue=similarities_by_group.columns[1],
                data=similarities_by_group,
                legend=False,
                palette=SUBPOCKET_COLORS,
            )
        except KeyError:
            sns.boxplot(
                x=similarities_by_group.columns[1],
                y=similarities_by_group.columns[0],
                data=similarities_by_group,
                color="dodgerblue",
            )
        plt.ylabel("Tanimoto similarity", fontsize=18)
        plt.xlabel(group_name, fontsize=18)
        plt.title(library_names[i])
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        i = i + 1

    plt.show()


def most_common_in_subset(
    fragment_library,
    clustered_fragments_subset,
    subpocket,
    top_x=50,
    num_subpockets=6,
):
    """
    Find the top_x most common fragments for the specified subpocket in the original,
    not deduplicated, fragment library and compare them with the clustered fragment library subset,
    getting the number of occurrences in the original fragment library, the cluster id from the
    clustered fragment library subset and the number of subpockets where
    the fragment is found in the subset.
    Draw fragments sorted by descending cluster size and subpocket count.

    Parameters
    ----------
    fragment_library : dict
        Deduplicated fragment library organized in subpockets.
    clustered_fragments_subset : pandas DataFrame
        Clustered fragments from a subset of the fragment_library
    subpocket : str
        Subpocket for which the top_x fragments are calculated.
    top_x : int
        Number of fragments with the highest occurrence that should be used for comparison.
        By default, top_x=50.
    num_subpockets : int
        Number of subpockets a fragment should at most be found in. Can be an integer value
        between 1 and 6.
        By default, num_subpockets=6.

    Returns
    -------
    PIL.PngImagePlugin.PngImageFile
        Image of fragments sorted by descending cluster size and subpocket count.

    """

    # if no fragments found that are in <=num_subpockets take one subpocket more (e.g. if no
    # fragment in 1 subpocket found take fragments in 2 subpockets)
    for i in range(num_subpockets, 7, 1):
        subset = clustered_fragments_subset[
            clustered_fragments_subset["subpocket_count"] <= i
        ]
        subset = subset[subset["subpockets"].astype(str).str.contains(subpocket)]
        if subset.empty:
            i = i + 1
        else:
            break
    # get the 50 most common fragments from the original fragment library which is
    # not deduplicted...
    most_common_frags = kfl_utils.get_most_common_fragments(
        fragment_library[subpocket],
        top_x=top_x,
    )
    # sort the most common fragments by the number of occourences of each fragment
    most_common_frags.sort_values(by=["fragment_count"], ascending=False, inplace=True)

    # compare the most common fragments from the original fragment library with the subset and keep
    # the ones that are in both
    most_common_in_subset = []
    cluster_ids = []
    fragment_counts = []
    most_common_subset_df = pd.DataFrame()
    subpocket_counts = []
    legend_lst = []
    for _, row_common in most_common_frags.iterrows():
        for _, row_subset in subset.iterrows():
            if Chem.MolFromSmiles(row_common["smiles"]).HasSubstructMatch(
                Chem.MolFromSmiles(row_subset["smiles"])
            ) and Chem.MolFromSmiles(row_subset["smiles"]).HasSubstructMatch(
                Chem.MolFromSmiles(row_common["smiles"])
            ):
                most_common_in_subset.append(row_subset["smiles"])
                fragment_counts.append(row_common["fragment_count"])
                cluster_ids.append(row_subset["cluster_id"])
                subpocket_counts.append(row_subset["subpocket_count"])

    most_common_subset_df["smiles"] = most_common_in_subset
    most_common_subset_df["cluster_id"] = cluster_ids
    most_common_subset_df["fragment_count"] = fragment_counts
    most_common_subset_df["subpocket_count"] = subpocket_counts

    # if there are fragments found in both subsets
    if not most_common_subset_df.empty:
        # sort the fragments by clster_id ascending and fragment_count descending
        most_common_subset_df.sort_values(
            by=["cluster_id", "fragment_count"], ascending=[True, False], inplace=True
        )
        # create the legend list with cluster_id, fragment_count and subpocket_count
        for _, row in most_common_subset_df.iterrows():
            legend_lst.append(
                f'{row["cluster_id"]} | {row["fragment_count"]} | {row["subpocket_count"]}'
            )
        # draw molecules
        print(
            f"Legend: subset cluster ID | fragment count inside %s in complete fragment library | fragment subpocket count in subset"
            % subpocket
        )  # noqa E501
        img = Draw.MolsToGridImage(
            [Chem.MolFromSmiles(smiles) for smiles in most_common_subset_df["smiles"]],
            legends=legend_lst,
            molsPerRow=10,
        )
        return img
    else:
        print("No fragment was found in both subsets")
