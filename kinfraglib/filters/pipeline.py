"""
Contains function to start the filters using the pipeline
"""
from kinfraglib import filters
import pandas as pd
from IPython.display import display
import warnings


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
    PATH_DATA_CUSTOM = general_parameters.get("custom_path")
    save_cols = []      # variable to store filter columns that are created during filtering
    if pains_parameters.get('pains_filter'):
        save_cols.append("bool_pains")
        print("Apply PAINS filter..")
        fragment_library, pains_df = filters.unwanted_substructures.get_pains(fragment_library)
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
        save_cols.append("bool_brenk")
        print("Apply Brenk filter..")
        DATA_BRENK = brenk_parameters.get("substructure_file")
        fragment_library, brenk_structs = filters.unwanted_substructures.get_brenk(
            fragment_library, DATA_BRENK
        )

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
        save_cols.append("bool_ro3")
        print("Apply Ro3 filter..")
        min_fulfilled_ro3 = ro3_parameters.get("min_fulfilled")
        cutoff_crit_ro3 = ro3_parameters.get("cutoff_crit")
        fragment_library = filters.drug_likeness.get_ro3_frags(
            fragment_library,
            min_fulfilled=min_fulfilled_ro3,
            cutoff_crit=cutoff_crit_ro3,
        )

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
        save_cols.append("bool_qed")
        save_cols.append("qed")
        print("Apply QED filter..")
        cutoff_qed = qed_parameters.get("cutoff_value")
        cutoff_crit_qed = qed_parameters.get("cutoff_crit")
        plot_stats_qed = qed_parameters.get("plot_stats")

        fragment_library = filters.drug_likeness.get_qed(
            fragment_library,
            cutoff_val=cutoff_qed,
            cutoff_crit=cutoff_crit_qed,
        )
        if general_parameters.get("show_stats"):
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
                fragment_library, "qed", "QED", plot_stats=plot_stats_qed, cutoff=cutoff_qed
            )

    if bb_parameters.get('bb_filter'):
        save_cols.append("bool_bb")
        print("Apply BB filter..")
        fragment_library = filters.synthesizability.check_building_blocks(
            fragment_library,
            bb_parameters.get("bb_file"),
        )
        if general_parameters.get("show_stats"):
            num_fragments_bb = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_bb", "enamine"
                    ),
                ],
                axis=1,
            )
            num_fragments_bb.append(num_fragments_bb.sum().rename('Total'))
            display(num_fragments_bb)
    if syba_parameters.get('syba_filter'):
        save_cols.append("bool_syba")
        save_cols.append("syba")
        print("Apply SYBA filter..")
        cutoff_syba = syba_parameters.get("cutoff_value")
        cutoff_crit_syba = syba_parameters.get("cutoff_crit")
        query_type = syba_parameters.get("query_type")
        plot_stats_syba = syba_parameters.get("plot_stats")
        fragment_library = filters.synthesizability.calc_syba(
            fragment_library,
            cutoff=cutoff_syba,
            cutoff_criteria=cutoff_crit_syba,
            query_type=query_type,
        )
    if general_parameters.get("show_stats"):
        num_fragments_syba = pd.concat(
            [
                filters.analysis.count_fragments(
                    fragment_library, "pre_filtered"
                ),
                filters.analysis.count_accepted_rejected(
                    fragment_library, "bool_syba", "syba"
                ),
            ],
            axis=1,
        )
        num_fragments_syba.append(num_fragments_syba.sum().rename('Total'))
        display(num_fragments_syba)
    if syba_parameters.get("do_plot"):
        filters.plots.make_hists(
            fragment_library,
            colname="syba",
            filtername="SYBA",
            plot_stats=plot_stats_syba,
            cutoff=cutoff_syba
        )
    # save filter results up to this point
    filters.retro.save_filter_results(fragment_library, save_cols, PATH_DATA_CUSTOM)
    fragment_library = filters.analysis.number_of_accepted(
        fragment_library, columns=[
            "bool_pains",
            "bool_brenk",
            "bool_ro3",
            "bool_qed",
            "bool_bb",
            "bool_syba",
        ],
        min_accepted=6,
    )

    for subpocket in fragment_library.keys():
        fragment_library[subpocket].drop(
            fragment_library[subpocket].loc[fragment_library[subpocket]['bool'] == 0].index,
            inplace=True,
        )
        fragment_library[subpocket] = fragment_library[subpocket].reset_index(drop=True)

    if retro_parameters.get('retro_filter'):
        print("Apply pairwise retrosynthesizability filter..")
        PATH_DATA_RETRO = retro_parameters.get("retro_path")
        valid_fragment_pairs, unique_pairs = filters.retro.get_valid_fragment_pairs(
            fragment_library
        )

        warnings.filterwarnings("ignore")
        fragment_library, mol_df, diff_df = filters.retro.get_pairwise_retrosynthesizability(
            unique_pairs,
            PATH_DATA_RETRO,
            valid_fragment_pairs,
            fragment_library,
        )

        if retro_parameters.get("plot_stats"):
            num_fragments_retro = pd.concat(
                [
                    filters.analysis.count_fragments(
                        fragment_library, "custom_filtered"
                    ),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_retro", "pairwise_retosynthesizability"
                    ),
                ],
                axis=1,
            )
            num_fragments_retro.append(num_fragments_retro.sum().rename("Total"))
            display(num_fragments_retro)

        if retro_parameters.get("do_plot"):
            filters.plots.make_retro_hists(fragment_library, "retro_count", cutoff=0)

        if retro_parameters.get("show_mols"):
            for subpocket in fragment_library.keys():
                plt1 = filters.plots.retro_routes_fragments(
                    fragment_library,
                    evaluate="none",
                    subpocket=subpocket,
                    molsPerRow=10,
                )
                plt1.show()
                plt2 = filters.plots.retro_routes_fragments(
                    fragment_library,
                    evaluate="max",
                    subpocket=subpocket,
                    molsPerRow=10,
                )
                plt2.show()

        saved_filter_results = pd.read_csv(PATH_DATA_CUSTOM / "custom_filter_results.csv")
        saved_filter_results.set_index(["subpocket", "smiles"])
        fragment_library_concat = pd.concat(fragment_library)
        retro_results_df = pd.DataFrame()
        retro_results_df["subpocket"] = fragment_library_concat["subpocket"]
        retro_results_df["smiles"] = fragment_library_concat["smiles"]
        retro_results_df["retro_count"] = fragment_library_concat["retro_count"]
        retro_results_df["bool_retro"] = fragment_library_concat["bool_retro"]
        retro_results_df.set_index(["subpocket", "smiles"])
        all_results_df = saved_filter_results.merge(
            retro_results_df,
            left_on=['subpocket', 'smiles'],
            right_on=['subpocket', 'smiles'],
            how="outer")
        all_results_df.to_csv(PATH_DATA_CUSTOM / "custom_filter_results.csv", index=False)

    print("Save custom filtered fragment library..")
    for subpocket in fragment_library.keys():
        fragment_library[subpocket].drop(
            fragment_library[subpocket].loc[fragment_library[subpocket]["bool_retro"] == 0].index,
            inplace=True,
        )
    filters.retro.save_fragment_library_to_sdfs(
        PATH_DATA_CUSTOM,
        fragment_library_concat,
    )
    fragment_library[subpocket] = fragment_library[subpocket].reset_index(drop=True)
    return {
        'fragment_library': fragment_library,
        'pains_df': pains_df,
        'brenk_df': brenk_structs,
        "mol_df": mol_df,
        "diff_df": diff_df,
    }
