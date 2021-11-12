"""
Contains plotting functionalities for different filters.
"""
import statistics
import matplotlib.pyplot as plt
from matplotlib import gridspec


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
                    fragment_library[keys[subpocket_num]][colname], facecolor="#04D8B2",
                    edgecolor="#808080"
                )
                ax.set_title(keys[((i * 4) + j)])
                if plot_stats:      # add statistics box (max, min, mean value per subpocket)
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="mean: " +    # noqa: W504
                        str(round(statistics.mean(fragment_library[keys[subpocket_num]][colname]))),     # noqa: E501
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="min: " +     # noqa: W504
                        str(round(min(fragment_library[keys[subpocket_num]][colname]))),
                    )
                    plt.plot(
                        [],
                        [],
                        " ",
                        label="max: " +  # noqa: W504
                        str(round(max(fragment_library[keys[subpocket_num]][colname]))),
                    )
                    plt.legend()
                if cutoff is not None:      # if a cutoff is given draw a red line
                    plt.axvline(x=cutoff, color="r", linestyle="-")
                if filtername is not None:
                    plt.xlabel(filtername)
                plt.ylabel("Number of fragments")       # set yaxis label
                subpocket_num = subpocket_num + 1       # go to next subpocket
    plt.suptitle(filtername)        # set filtername as title over the plots
    plt.show()      # show the plots, needed when called in function other not in notebook
