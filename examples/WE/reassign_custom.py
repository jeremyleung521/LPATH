import numpy

def reassign_custom(data, pathways, dictionary, assign_file=None):
    """
    Reclassify/assign frames into different states. This is highly
    specific to the system. If w_assign's definition is sufficient,
    you can proceed with what's made in the previous step
    using ``reassign_identity``.

    In this example, the dictionary maps state idx to its corresponding ``state_string``.
    We suggest using alphabets as states.

    Parameters
    ----------
    data : list
        An array with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``. The last entry in
        the dictionary should be the "unknown" state.

    assign_file : str, default : None
        A string pointing to the ``assign.h5`` file. Needed as a parameter for all functions,
        but is ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a ``state string`` (character).
        The last entry in the dictionary should be the "unknown" state.

    """
    # reassign states to be the cluster IDs
    for idx, val in enumerate(data):  # Loop through each set of successful pathways
        val_arr = numpy.asarray(val)
        for idx2, val2 in enumerate(val_arr):  # Loop through each frame of the pathway
                val2[2] = int(val2[-3])  # Renumber state_id the with the aux dataset
                pathways[idx, idx2] = val2

    # Generating a dictionary mapping each state
    dictionary = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6", 7: "7"}

    return dictionary
