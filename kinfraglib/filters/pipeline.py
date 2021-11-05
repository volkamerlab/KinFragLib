"""
Contains function to start the filters using the pipeline
"""
from kinfraglib import filters
import pandas as pd
from IPython.display import display


def start_pipeline(
    fragment_library,
    pains_parameters,
    brenk_parameters,
    ro3_parameters,
    qed_parameters,
    bb_parameters,
    syba_parameters,
    retro_parameters,
    general_parameters,
):
    if pains_parameters.get('pains_filter'):
        print("Apply PAINS filter.")
        pains_dict = filters.pains.get_pains(fragment_library)
        fragment_library = pains_dict["fragment_library"]
        pains_df = pains_dict["pains"]
        if general_parameters.get("show_stats"):
            num_fragments_pains = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_pains", "pains"
                    ),
                ],
                axis=1,
            )
            num_fragments_pains.append(num_fragments_pains.sum().rename('Total'))

            display(num_fragments_pains)

    if brenk_parameters.get("brenk_filter"):
        print("Apply Brenk filter.")
        DATA_BRENK = brenk_parameters.get("substructure_file")
        brenk_dict = filters.unwanted_substructures.get_brenk(
            fragment_library, DATA_BRENK
        )
        fragment_library = brenk_dict["fragment_library"]
        brenk_structs = brenk_dict["brenk"]

        if general_parameters.get("show_stats"):
            num_fragments_brenk = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_brenk", "brenk"
                    ),
                ],
                axis=1,
            )
            num_fragments_brenk.append(num_fragments_brenk.sum().rename('Total'))
            display(num_fragments_brenk)

    if ro3_parameters.get('ro3_filter'):
        print("Apply Ro3 filter.")
        fragment_library = filters.ruleofthree.get_ro3_frags(fragment_library)

        if general_parameters.get("show_stats"):
            num_fragments_ro3 = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_ro3", "ro3"
                    ),
                ],
                axis=1,
            )
            num_fragments_ro3.append(num_fragments_ro3.sum().rename('Total'))

            display(num_fragments_ro3)

    if qed_parameters.get('qed_filter'):
        print("Apply QED filter.")
        cutoff_qed = qed_parameters.get("cutoff_value")
        cutoff_crit_qed = qed_parameters.get("cutoff_crit")

        if general_parameters.get("show_stats"):
            fragment_library = filters.qed.get_qed(
                fragment_library,
                cutoff_val=cutoff_qed,
                cutoff_crit=cutoff_crit_qed,
            )
            num_fragments_qed = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_qed", "qed"
                    ),
                ],
                axis=1,
            )
            num_fragments_qed.append(num_fragments_qed.sum().rename('Total'))
            display(num_fragments_qed)

        if qed_parameters.get("do_plot"):
            filters.plots.make_hists(
                fragment_library, "qed", "QED", plot_stats=True, cutoff=cutoff_qed
            )

    if bb_parameters.get('bb_filter'):
        print("Apply BB filter.")
    if syba_parameters.get('syba_filter'):
        print("Apply SYBA filter.")
    if qed_parameters.get('retro_filter'):
        print("Apply pairwise retrosynthesizability filter.")

    return {
        'fragment_library': fragment_library,
        'pains_df': pains_df,
        'brenk_df': brenk_structs,
    }
