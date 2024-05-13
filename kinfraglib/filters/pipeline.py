"""
Contains function to start the filters using the pipeline
"""
from kinfraglib import filters
import pandas as pd
from IPython.display import display
import warnings
import os
import datetime
import pathlib
import json     # noqa F401
import copy     # noqa F401


def start_pipeline(
    fragment_library,
    pains_parameters,
    brenk_parameters,
    ro3_parameters,
    qed_parameters,
    bb_parameters,
    syba_parameters,
    retro_parameters,
    global_parameters,
):
    """
    Starting the custom filters' pipeline with the defined parameters

    Parameters
    ----------
    fragment_library : dict
        fragments organized in subpockets including all information

    pains_parameters : dict
        containing the following parameters: 'pains_filter' (boolean)

    brenk_parameters : dict
        containing the following parameters: 'brenk_filter'(boolean),
        'substructure_file_path'(Path)

    ro3_parameters : dict
        containing the following parameters: 'ro3_filter'(boolean), 'num_fulfilled'(int),
        `cutoff_crit`(str)

    qed_parameters : dict
        containing the following parameters: 'qed_filter'(boolean), 'cutoff_value'(float),
        'cutoff_crit'(str), 'do_plot'(boolean), 'plot_stats'(boolean)

    bb_parameters : dict
        containing the following parameters: 'bb_filter'(boolean), 'bb_file'(Path)

    syba_parameters : dict
        containing the following parameters: 'syba_filter'(boolean), 'cutoff_value'(int),
        'cutoff_crit'(str), 'query_type'(str), 'do_plot'(boolean), 'plot_stats(boolean)

    retro_parameters : dict
        containing the following parameters: 'retro_filter'(boolean), 'cutoff_value'(int),
        'cutoff_crit'(str), 'retro_path'(Path), 'do_plot'(boolean), 'show_mols'(boolean),
        'plot_stats'(boolean)

    global_parameters: dict
        containing the following parameters: 'show_stats'(boolean), 'custom_path'(Path),
        'num_passing'(int)

    Returns
    dict
        filtered fragment library
        data created during filtering
    -------

    """
    # check if number of filters to be passed is not bigger than actual number of filters applied
    num_filters = sum(
        [
            pains_parameters.get("pains_filter"),
            brenk_parameters.get("brenk_filter"),
            ro3_parameters.get("ro3_filter"),
            qed_parameters.get("qed_filter"),
            bb_parameters.get("bb_filter"),
            syba_parameters.get("syba_filter"),
        ]
    )
    num_passing = global_parameters.get("num_passing")
    if retro_parameters.get("retro_filter"):
        if num_filters < num_passing:
            print_str = (
                "Only %s filters are activated before applying pairwise retrosynthesizability. \n"
                "Setting `num_passing` to %s." % (num_filters, num_filters)
            )
            print(print_str)
            num_passing = num_filters
    dir = os.path.join(
        global_parameters.get("custom_path"),
        datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'),
    )

    if not os.path.exists(dir):
        os.makedirs(dir)

    PATH_DATA_CUSTOM = pathlib.Path(dir)    # define Path to custom data
    
    # to prevent showing personal paths in the notebooks
    dir_print_version = str(dir).split("../../")[1] if len(str(dir).split("../../")) == 2 else dir

    print(
        "Your custom kinfraglib, the chosen parameters log file and the filtering results will be stored in " 
        + dir_print_version   # noqa E501
    )
    # save chosen parameters in created timestamp dir to save log file and created library
    param_list = [
        "global_parameters",
        "pains_parameters",
        "brenk_parameters",
        "ro3_parameters",
        "qed_parameters",
        "bb_parameters",
        "syba_parameters",
        "retro_parameters",
    ]

    for curdict in param_list:
        cur_dict = eval("copy.deepcopy(" + curdict + ")")
        for key in eval(curdict + ".keys()"):
            if eval("isinstance(cur_dict[\"" + key + "\"], pathlib.Path)"):
                new_path = {key: str(cur_dict[key])}  # noqa F841
                cur_dict.update(new_path)
                cur_dict[key] = str(cur_dict[key])
        with open(PATH_DATA_CUSTOM / "custom_filtering_parameters.log", 'a+') as fp:
            fp.write(curdict + ": ")

        eval("json.dump(cur_dict, open(\"" + str(PATH_DATA_CUSTOM) + "/custom_filtering_parameters.log\", \"a+\"))")     # noqa E501
        with open(PATH_DATA_CUSTOM / "custom_filtering_parameters.log", 'a+') as fp:
            fp.write("\n")

    save_cols = []  # variable to store filter columns that are created during filtering

    # if pains_filter is activated, apply pains filter with the given parameters
    if pains_parameters.get("pains_filter"):
        save_cols.append("bool_pains")
        print("Apply PAINS filter..")
        fragment_library, pains_df = filters.unwanted_substructures.get_pains(
            fragment_library
        )
        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_pains = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_pains", "pains"
                    ),
                ],
                axis=1,
            )
            num_fragments_pains = pd.concat(
                [num_fragments_pains, num_fragments_pains.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_pains)
    # if brenk_filter is activated, apply brenk filter with the given parameters
    if brenk_parameters.get("brenk_filter"):
        save_cols.append("bool_brenk")
        print("Apply Brenk filter..")
        DATA_BRENK = brenk_parameters.get("substructure_file_path")
        fragment_library, brenk_structs = filters.unwanted_substructures.get_brenk(
            fragment_library, DATA_BRENK
        )

        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_brenk = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_brenk", "brenk"
                    ),
                ],
                axis=1,
            )
            num_fragments_brenk = pd.concat(
                [num_fragments_brenk, num_fragments_brenk.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_brenk)

    # if ro3_filter is activated, apply ro3 filter with the given parameters
    if ro3_parameters.get("ro3_filter"):
        save_cols.append("bool_ro3")
        print("Apply Ro3 filter..")
        num_fulfilled_ro3 = ro3_parameters.get("num_fulfilled")
        cutoff_crit_ro3 = ro3_parameters.get("cutoff_crit")
        fragment_library = filters.drug_likeness.get_ro3_frags(
            fragment_library,
            min_fulfilled=num_fulfilled_ro3,
            cutoff_crit=cutoff_crit_ro3,
        )

        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_ro3 = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_ro3", "ro3"
                    ),
                ],
                axis=1,
            )
            num_fragments_ro3 = pd.concat(
                [num_fragments_ro3, num_fragments_ro3.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_ro3)

    # if qed_filter is activated, apply qed filter with the given parameters
    if qed_parameters.get("qed_filter"):
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

        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_qed = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_qed", "qed"
                    ),
                ],
                axis=1,
            )
            num_fragments_qed = pd.concat(
                [num_fragments_qed, num_fragments_qed.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_qed)

        if qed_parameters.get("do_plot"):
            filters.plots.make_hists(
                fragment_library,
                "qed",
                "QED",
                plot_stats=plot_stats_qed,
                cutoff=cutoff_qed,
            )

    # if bb_filter is activated, apply bb filter with the given parameters
    if bb_parameters.get("bb_filter"):
        save_cols.append("bool_bb")
        print("Apply BB filter..")
        fragment_library = filters.synthesizability.check_building_blocks(
            fragment_library,
            bb_parameters.get("bb_file"),
        )

        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_bb = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_bb", "enamine"
                    ),
                ],
                axis=1,
            )
            num_fragments_bb = pd.concat(
                [num_fragments_bb, num_fragments_bb.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_bb)

    # if syba_filter is activated, apply syba filter with the given parameters
    if syba_parameters.get("syba_filter"):
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

        # if user wants to see statistics, plot number of fragments accepted/rejected
        if global_parameters.get("show_stats"):
            num_fragments_syba = pd.concat(
                [
                    filters.analysis.count_fragments(fragment_library, "pre_filtered"),
                    filters.analysis.count_accepted_rejected(
                        fragment_library, "bool_syba", "syba"
                    ),
                ],
                axis=1,
            )
            num_fragments_syba = pd.concat(
                [num_fragments_syba, num_fragments_syba.sum().rename("Total").to_frame().T]
            )
            display(num_fragments_syba)
        if syba_parameters.get("do_plot"):
            filters.plots.make_hists(
                fragment_library,
                colname="syba",
                filtername="SYBA",
                plot_stats=plot_stats_syba,
                cutoff=cutoff_syba,
            )
    # check which filters were applied
    nfrags = filters.analysis.count_fragments(fragment_library, "pre_filtered")
    # save filter results up to this point
    filters.retro.save_filter_results(fragment_library, save_cols, PATH_DATA_CUSTOM)
    filter_cols = []
    if pains_parameters.get("pains_filter"):
        filter_cols.append("bool_pains")
    if brenk_parameters.get("brenk_filter"):
        filter_cols.append("bool_brenk")
    if ro3_parameters.get("ro3_filter"):
        filter_cols.append("bool_ro3")
    if qed_parameters.get("qed_filter"):
        filter_cols.append("bool_qed")
    if bb_parameters.get("bb_filter"):
        filter_cols.append("bool_bb")
    if syba_parameters.get("syba_filter"):
        filter_cols.append("bool_syba")
    # get boolean column defining if a fragment passed at least the defined number of filters
    fragment_library = filters.analysis.number_of_accepted(
        fragment_library,
        columns=pd.Series(filter_cols),
        min_accepted=num_passing,
    )

    # save only fragments passing the defined min. number of filters
    for subpocket in fragment_library.keys():
        fragment_library[subpocket].drop(
            fragment_library[subpocket]
            .loc[fragment_library[subpocket]["bool"] == 0]
            .index,
            inplace=True,
        )
        fragment_library[subpocket] = fragment_library[subpocket].reset_index(drop=True)

    # create a dataframe counting how many fragments are used for pairwise retro filter
    frags_for_retro = pd.DataFrame(
        pd.concat(
            [
                nfrags,
                filters.analysis.count_fragments(
                    fragment_library,
                    "number of fragments used for pairwise retrosynthesizability",
                ),
            ],
            axis=1,
        )
    )
    frags_for_retro = pd.concat(
        [frags_for_retro, frags_for_retro.sum().rename("Total").to_frame().T]
    )
    display(frags_for_retro)

    if retro_parameters.get("retro_filter"):
        print("Apply pairwise retrosynthesizability filter..")
        PATH_DATA_RETRO = retro_parameters.get("retro_path")
        valid_fragment_pairs, unique_pairs = filters.retro.get_valid_fragment_pairs(
            fragment_library
        )

        warnings.filterwarnings("ignore")
        (
            fragment_library,
            mol_df,
            diff_df,
        ) = filters.retro.get_pairwise_retrosynthesizability(
            unique_pairs,
            PATH_DATA_RETRO,
            valid_fragment_pairs,
            fragment_library,
            cutoff_value=retro_parameters.get("cutoff_value"),
            cutoff_crit=retro_parameters.get("cutoff_crit"),
        )

        # if user wants to see statistics, plot number of fragments accepted/rejected
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
            num_fragments_retro = num_fragments_retro.append(
                num_fragments_retro.sum().rename("Total")
            )
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
                display(plt1)
                plt2 = filters.plots.retro_routes_fragments(
                    fragment_library,
                    evaluate="max",
                    subpocket=subpocket,
                    molsPerRow=10,
                )
                display(plt2)

        # load results from filtering steps and add the pairwise retrosynthesizability results
        saved_filter_results = pd.read_csv(
            PATH_DATA_CUSTOM / "custom_filter_results.csv"
        )
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
            left_on=["subpocket", "smiles"],
            right_on=["subpocket", "smiles"],
            how="outer",
        )
        all_results_df.to_csv(
            PATH_DATA_CUSTOM / "custom_filter_results.csv", index=False
        )

    # save the filtered fragment library
    print("Save custom filtered fragment library to %s" % str(PATH_DATA_CUSTOM))
    for subpocket in fragment_library.keys():
        fragment_library[subpocket].drop(
            fragment_library[subpocket]
            .loc[fragment_library[subpocket]["bool_retro"] == 0]
            .index,
            inplace=True,
        )
    # remove fragments not passing the retro-filter
    fragment_library_concat = fragment_library_concat[fragment_library_concat["bool_retro"] == 1]
    filters.retro.save_fragment_library_to_sdfs(
        PATH_DATA_CUSTOM,
        fragment_library_concat,
    )
    fragment_library[subpocket] = fragment_library[subpocket].reset_index(drop=True)
    return {
        "fragment_library": fragment_library,
        "pains_df": pains_df,
        "brenk_df": brenk_structs,
        "mol_df": mol_df,
        "diff_df": diff_df,
    }
