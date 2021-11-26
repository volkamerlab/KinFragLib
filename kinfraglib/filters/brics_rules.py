# copied from
#  https://github.com/volkamerlab/KinaseFocusedFragmentLibrary/blob/b7e684c26f75efffc2a9ba2383c9027cdd4c29a3/kinase_focused_fragment_library/recombination/classes_meta.py**  # noqa: E501

# BRICS rules as defined by RDKit
brics_rules = [
    {"1", "3"},
    {"1", "5"},
    {"1", "10"},
    # L2 definition is incorporated into L5
    {"3", "4"},
    {"3", "13"},
    {"3", "14"},
    {"3", "15"},
    {"3", "16"},
    {"4", "5"},
    {"4", "11"},
    {"5", "12"},
    {"5", "13"},
    {"5", "14"},
    {"5", "15"},
    {"5", "16"},
    {"6", "13"},
    {"6", "14"},
    {"6", "15"},
    {"6", "16"},
    {"7", "7"},
    {"8", "9"},
    {"8", "10"},
    {"8", "13"},
    {"8", "14"},
    {"8", "15"},
    {"8", "16"},
    {"9", "13"},  # not in original paper
    {"9", "14"},  # not in original paper
    {"9", "15"},
    {"9", "16"},
    {"10", "13"},
    {"10", "14"},
    {"10", "15"},
    {"10", "16"},
    {"11", "13"},
    {"11", "14"},
    {"11", "15"},
    {"11", "16"},
    {"13", "14"},
    {"13", "15"},
    {"13", "16"},
    {"14", "14"},  # not in original paper
    {"14", "15"},
    {"14", "16"},
    {"15", "16"},
    {"16", "16"},  # not in original paper
]


def is_brics_bond(env_1, env_2):

    """
    Checks if two given BRICS environment types are allowed to be connected according to the BRICS
    algorithm
    Parameters
    ----------
    env_1, env_2: str
        BRICS environment types
    Returns
    -------
    True if the given environments are allowed to be connected
    False otherwise
    """

    if {env_1, env_2} in brics_rules:
        return True

    else:
        return False
