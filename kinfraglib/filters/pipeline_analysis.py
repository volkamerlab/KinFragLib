"""
Contains functions for analysing the libraries (prefiltered vs. reduced vs. custom)
"""
from kinfraglib import utils as kfl_utils
from rdkit.Chem import Draw
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


SUBPOCKET_COLORS = {
    "AP": "purple",
    "FP": "forestgreen",
    "SE": "c",
    "GA": "tab:orange",
    "B1": "tab:blue",
    "B2": "darkslateblue",
    "X": "grey",
}


def get_clustered_most_common_fragments(fragment_library, subpocket):
    # function copied from https://github.com/volkamerlab/KinFragLib/blob/master/notebooks/2_3_fragment_analysis_most_common_fragments.ipynb  # noqa: E501
    """
    Get top X (50 per default) most common fragments (if multiple fragments have the same count
    but some make it into the top X and some not, include the latter also).
    Cluster the most common fragments by fingerprint similarity.

    Parameters
    ----------
    fragment_library : dict of pandas.DataFrame
        Fragment details (values), i.e. SMILES, kinase groups, and fragment RDKit molecules,
        for each subpocket (key).
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.

    Returns
    -------
    pandas.DataFrame
        Most common fragments (ID, SMILES, ROMol, cluster ID, fragment count).
    """

    # Get top X most common fragments
    if subpocket == "all":
        most_common_fragments = kfl_utils.get_most_common_fragments(
            pd.concat(fragment_library),
            top_x=50,
        )
    else:
        most_common_fragments = kfl_utils.get_most_common_fragments(
            fragment_library[subpocket],
            top_x=50,
        )

    # Cluster fingerprints
    clusters = kfl_utils.cluster_molecules(most_common_fragments.ROMol, cutoff=0.6)

    # Link fragments to cluster ID
    most_common_fragments = most_common_fragments.merge(
        clusters,
        on='molecule_id'
    )

    most_common_fragments.sort_values(
        ['cluster_id', 'fragment_count'],
        ascending=[True, False],
        inplace=True
    )

    most_common_fragments.reset_index(inplace=True, drop=True)

    return most_common_fragments


def plot_cluster_sizes(most_common_fragments, subpocket):
    # function copied from https://github.com/volkamerlab/KinFragLib/blob/master/notebooks/2_3_fragment_analysis_most_common_fragments.ipynb  # noqa: E501
    """
    Plot cluster sizes (cluster ID vs. cluster size).

    Parameters
    ----------
    most_common_fragments : pandas.DataFrame
        Most common fragments (ID, SMILES, ROMol, cluster ID, fragment count).
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.
    """

    cluster_sizes = most_common_fragments.groupby('cluster_id').size()
    cluster_sizes.name = 'cluster_size'
    cluster_sizes.plot(kind='bar', title=f'Cluster sizes for subpocket {subpocket}')


def draw_clusters(most_common_fragments, subpocket, output_path=None):
    # function copied from https://github.com/volkamerlab/KinFragLib/blob/master/notebooks/2_3_fragment_analysis_most_common_fragments.ipynb  # noqa: E501
    """
    Draw fragments sorted by descending cluster size and fragment count.

    Parameters
    ----------
    most_common_fragments : pandas.DataFrame
        Most common fragments (ID, SMILES, ROMol, cluster ID, fragment count).
    subpocket : str
        Subpocket name, i.e. AP, SE, FP, GA, B1, or B2.
    output_path : pathlib.Path
        Path to output folder.

    Returns
    -------
    PIL.PngImagePlugin.PngImageFile
        Image of fragments sorted by descending cluster size.
    """

    img = Draw.MolsToGridImage(
        list(most_common_fragments.ROMol),
        legends=[
            f'{row.cluster_id} | {row.fragment_count}'
            for index, row
            in most_common_fragments.iterrows()
        ],
        molsPerRow=7,
        maxMols=100,
        subImgSize=(170, 170),
        useSVG=True
    )

    print('Legend: cluster ID | fragment count')

    if output_path is not None:

        # Get SVG data
        molsvg = img.data

        # Set font size
        molsvg = molsvg.replace('12px', '24px')

        # Save altered SVG data to file
        with open(
            Path(output_path) / f'clustered_most_common_fragments_{subpocket.lower()}.svg',
            'w',
        ) as f:
            f.write(molsvg)

    return img


def plot_fragment_descriptors(descriptors):
    """
    Plot fragment descriptors.
    """
    # copied from utils without saving img

    plt.figure(figsize=(25, 6))

    for i, descriptor_name in enumerate(descriptors.columns[3:]):

        plt.subplot(1, 4, i + 1)
        sns.boxplot(
            x="subpocket",
            y=descriptor_name,
            data=descriptors,
            palette=SUBPOCKET_COLORS,
            medianprops={"linewidth": 3, "linestyle": "-"},
        )
        plt.ylabel(descriptor_name, fontsize=16)
        plt.xlabel("Subpocket", fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
    return plt


def plot_fragment_similarity(similarities_by_groups, group_name):
    """
    Plot fragment similarity by category, such as subpocket or kinase group.
    """
    # copied from kinfraglib/utils.py without saving img

    plt.figure(figsize=(20, 5))
    num_plots = len(similarities_by_groups)
    i = 0
    for similarities_by_group in similarities_by_groups:
        plt.subplot(1, num_plots, i + 1)
        try:
            sns.boxplot(
                x=similarities_by_group.columns[1],
                y=similarities_by_group.columns[0],
                data=similarities_by_group,
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
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        i = i + 1

    plt.show()
