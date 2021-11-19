"""
Contains plotting functionalities for different filters.
"""
import statistics
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
from rdkit.Chem import Draw, MACCSkeys
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
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
        name of the filter used as title creating the values plotted. By default None, meaning
        no title will be displayed.

    plot_stats : boolean
        defining if a box with min, max and mean value will be displayed in each plot. By default
        plot_stats=True.

    cutoff : int or float
        cutoff value for drawing a cutoff line to the plots. By default cutoff=None, meaning no
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
            if (i * 4) + j < num_plots:
                ax = plt.subplot(gs[i, j])
                ax.hist(
                    fragment_library[keys[subpocket_num]][colname],
                    facecolor="#04D8B2",
                    edgecolor="#808080",
                )
                ax.set_title(keys[((i * 4) + j)])
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
    # get even number if number of plots not even
    num_plots = round(len(fragment_library.keys()) + 0.5)
    plt.figure(figsize=(25, 29))
    gs = gridspec.GridSpec(int(num_plots / 2), int(num_plots / 2))
    keys = list(fragment_library.keys())
    subpocket_num = 0
    for i in range(0, 2):
        for j in range(0, int((num_plots) / 2)):
            if (i * 4) + j <= num_plots:
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
        )
        return img


def create_tsne_plots(fragment_library):
    """
    Creates t-SNE plots comparing
    a) pre-filtered and reduced fragment library
    b) pre-filtered and custom filtered fragment library
    c) pre-filtered, reduced and custom fragment library

    and prints number of fragments in the subsets.
    ----------
    fragment_library : dict
        fragment library organized in subpockets containing boolean columuns `bool_reduced`and
        `bool_custom`defining if the fragments are part of the subsets

    """
    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    fragment_library_concat["maccs"] = fragment_library_concat.ROMol.apply(MACCSkeys.GenMACCSKeys)

    pca = PCA(n_components=50)
    crds = pca.fit_transform(list(fragment_library_concat["maccs"]))

    crds_embedded = TSNE(n_components=2, init='pca', learning_rate='auto').fit_transform(crds)

    tsne_df = pd.DataFrame(crds_embedded, columns=["X", "Y"])
    # add bool column from filtering steps here
    tsne_df['reduced'] = fragment_library_concat["bool_reduced"]
    tsne_df['custom'] = fragment_library_concat["bool_custom"]
    # create column defining if fragment is
    # *excluded in both subsets (0)
    # *excluded in reduced (1)
    # *excluded in custom (2)
    # *accepted in both subsets (3)
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
    fig = plt.figure(figsize=(13, 10))
    plt.subplot(2, 2, 1)
    fig.add_subplot(2, 2, 1)
    sns.scatterplot(
        data=tsne_df.query("reduced == 0"),
        x="X",
        y="Y",
        color='lightcoral',
        alpha=0.5
    ).set_title("pre_filtered vs. reduced")
    sns.scatterplot(
        data=tsne_df.query("reduced == 1"),
        x="X",
        y="Y",
        color='green',
        alpha=0.5)
    plt.axis('off')

    plt.subplot(2, 2, 2)
    fig.add_subplot(2, 2, 2)
    sns.scatterplot(
        data=tsne_df.query("custom == 0"),
        x="X",
        y="Y",
        color='lightcoral',
        alpha=0.5
    ).set_title("pre-filtered vs. custom")
    sns.scatterplot(
        data=tsne_df.query("custom == 1"),
        x="X",
        y="Y",
        color='green',
        alpha=0.5
    )
    plt.axis('off')

    plt.subplot(2, 2, 3)
    fig.add_subplot(2, 2, 3)
    sns.scatterplot(
        data=tsne_df.query("compare == 0"),
        x="X",
        y="Y",
        color='lightcoral',
        alpha=0.5
    ).set_title("pre-filtered vs. reduced vs. custom")
    sns.scatterplot(
        data=tsne_df.query("compare == 1"),
        x="X",
        y="Y",
        color='orange',
        alpha=0.5
    )
    sns.scatterplot(
        data=tsne_df.query("compare == 2"),
        x="X",
        y="Y",
        color='lightblue',
        alpha=0.5
    )
    sns.scatterplot(
        data=tsne_df.query("compare == 3"),
        x="X",
        y="Y",
        color='green',
        alpha=0.5
    )

    plt.axis('off')
    plt.show()
    num_lists = (len(tsne_df["compare"]), num0, num1, num2, num3)
    print("""%s Pre-filtered fragments.
        Number of fragments excluded in both datasets: %s
        Number of fragments excluded in the reduced dataset but included in the custom dataset: %s
        Number of fragments excluded in the custom dataset but included in the reduced dataset: %s
        Number of fragments in both datasets: %s """ % (num_lists))


def create_tsne_plots_filters(fragment_library, saved_filter_results):
    """
    Creates t-SNE plots with accepted (green) and rejected (red) fragments for each filtering step.

    ----------
    fragment_library : dict
        fragment library organized in subpockets containing boolean columuns
    saved_filter_results : dataframe
        loaded file with saved filter results

    """
    fragment_library_concat = pd.concat(fragment_library).reset_index(drop=True)
    fragment_library_concat["maccs"] = fragment_library_concat.ROMol.apply(MACCSkeys.GenMACCSKeys)

    pca = PCA(n_components=50)
    crds = pca.fit_transform(list(fragment_library_concat["maccs"]))

    crds_embedded = TSNE(n_components=2, init='pca', learning_rate='auto').fit_transform(crds)

    tsne_df = pd.DataFrame(crds_embedded, columns=["X", "Y"])
    # add bool column from filter steps
    filters = []
    if "bool_pains" in saved_filter_results.columns:
        tsne_df['pains'] = saved_filter_results["bool_pains"]
        filters.append("pains")
    if "bool_brenk" in saved_filter_results.columns:
        tsne_df['brenk'] = saved_filter_results["bool_brenk"]
        filters.append("brenk")
    if "bool_ro3" in saved_filter_results.columns:
        tsne_df['ro3'] = saved_filter_results["bool_ro3"]
        filters.append("ro3")
    if "bool_qed" in saved_filter_results.columns:
        tsne_df["qed"] = saved_filter_results["bool_qed"]
        filters.append("qed")
    if "bool_bb" in saved_filter_results.columns:
        tsne_df['bb'] = saved_filter_results["bool_bb"]
        filters.append("bb")
    if "bool_syba" in saved_filter_results.columns:
        tsne_df['syba'] = saved_filter_results["bool_syba"]
        filters.append("syba")
    if "bool_retro" in saved_filter_results.columns:
        tsne_df['retro'] = saved_filter_results["bool_retro"]
        filters.append("retro")

    # create the plots for all filters
    fig = plt.figure(figsize=(15, 15))
    i = 0
    for filter in filters:
        i = i + 1
        plt.subplot(4, 2, i)
        fig.add_subplot(4, 2, i)
        sns.scatterplot(
            data=tsne_df.query("%s == 1" % filter),
            x="X",
            y="Y",
            color='green',
            alpha=0.5,
        ).set_title(filter)
        sns.scatterplot(data=tsne_df.query("%s == 0" % filter), x="X", y="Y", color='lightcoral')
        plt.axis('off')


def connection_frequencies(fragment_library, fragment_library_reduced, fragment_library_custom):
    """
    Calculates the connection frequencies between the subpockets for all three subsets and
    creates a plot.

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
    connections_by_fragment = kfl_utils.get_connections_by_fragment(fragment_library_concat)
    connections_by_ligand = connections_by_fragment.groupby(
        ['kinase', 'complex_pdb', 'ligand_pdb']
    )['connections_name'].sum()
    connections_by_ligand_count = connections_by_ligand.apply(lambda x: Counter(x))
    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count = pd.Series(
        Counter(connections_by_ligand_count.apply(list).sum())
    )
    connections_across_ligands_count.name = 'count_pre-filtered'

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency = connections_across_ligands_count.apply(
        lambda x: round((x / connections_by_ligand_count.shape[0] * 100), 1)
    )
    connections_across_ligands_frequency.name = 'frequency_pre-filtered'

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
        ['kinase', 'complex_pdb', 'ligand_pdb']
    )['connections_name'].sum()
    connections_by_ligand_count_reduced = connections_by_ligand_reduced.apply(lambda x: Counter(x))

    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count_reduced = pd.Series(
        Counter(connections_by_ligand_count_reduced.apply(list).sum())
    )
    connections_across_ligands_count_reduced.name = 'count_reduced'

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency_reduced = connections_across_ligands_count_reduced.apply(
        lambda x: round((x / connections_by_ligand_count_reduced.shape[0] * 100), 1)
    )
    connections_across_ligands_frequency_reduced.name = 'frequency_reduced'

    # Concatenate count and frequency data to DataFrame
    connections_across_ligands_reduced = pd.concat(
        [connections_across_ligands_count_reduced, connections_across_ligands_frequency_reduced],
        axis=1,
    )

    # fragment library custom filtered
    fragment_library_custom_concat = pd.concat(fragment_library_custom)
    connections_by_fragment_custom = kfl_utils.get_connections_by_fragment(
        fragment_library_custom_concat
    )

    connections_by_ligand_custom = connections_by_fragment_custom.groupby(
        ['kinase', 'complex_pdb', 'ligand_pdb']
    )['connections_name'].sum()
    connections_by_ligand_count_custom = connections_by_ligand_custom.apply(lambda x: Counter(x))

    # Get connection count across ligands (count each connection per ligand only once)
    connections_across_ligands_count_custom = pd.Series(
        Counter(connections_by_ligand_count_custom.apply(list).sum())
    )
    connections_across_ligands_count_custom.name = 'count_custom-filtered'

    # Get connection frequency (100% = all ligands)
    connections_across_ligands_frequency_custom = connections_across_ligands_count_custom.apply(
        lambda x: round((x / connections_by_ligand_count_custom.shape[0] * 100), 1)
    )
    connections_across_ligands_frequency_custom.name = 'frequency_custom_filtered'

    # Concatenate count and frequency data to DataFrame
    connections_across_ligands_custom = pd.concat(
        [connections_across_ligands_count_custom, connections_across_ligands_frequency_custom],
        axis=1,
    )

    frequencies = pd.concat(
        [
            connections_across_ligands["frequency_pre-filtered"],
            connections_across_ligands_reduced["frequency_reduced"],
            connections_across_ligands_custom["frequency_custom_filtered"],
        ],
        axis=1,
    )

    ax = frequencies.plot.bar()
    fig = ax.get_figure()

    fig.set_figheight(5)
    fig.set_figwidth(13)

    ax.set_xlabel("Subpocket")
    ax.set_ylabel("Connection Frequency")
    ax.set_title("Connection Frequencies of the different subsets")

    fig.show()
    res = pd.concat(
        [
            connections_across_ligands,
            connections_across_ligands_reduced,
            connections_across_ligands_custom
        ],
        axis=1,
    )
    return res


def num_frags_development(filter_res):
    frag_keys = filter_res.keys()
    frag_keys.to_list()
    bool_keys = [x for x in frag_keys if "bool" in x]
    update_results = pd.DataFrame()
    update_results["pre-filtered"] = filter_res.reset_index().groupby(
        "subpocket", sort=False
    ).size()
    for bool_key in bool_keys:
        filter_res = filter_res.loc[filter_res[bool_key] == 1]
        # filter_res = filter_res[filter_res[bool_key].notnull()]
        update_results[bool_key] = filter_res.reset_index().groupby("subpocket", sort=False).size()

    ax = update_results.plot.bar()
    fig = ax.get_figure()

    fig.set_figheight(5)
    fig.set_figwidth(13)

    ax.set_xlabel("Subpocket")
    ax.set_ylabel("Number of fragments")
    ax.set_title("Development of the number of fragments per subpocket after each filter step")

    fig.show()

    return update_results
