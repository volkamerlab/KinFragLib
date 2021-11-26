"""
Contains plotting functionalities for different filters.
"""
import statistics
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
from rdkit.Chem import Draw


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
