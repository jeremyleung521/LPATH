import numpy

def reassign_custom(data, pathways, dictionary, assign_file=None):

    # reassign states to be the cluster IDs
    for idx, val in enumerate(data):
        val_arr = numpy.asarray(val)
        for idx2, val2 in enumerate(val_arr):
                val2[2] = int(val2[3])
                pathways[idx, idx2] = val2

    # Generating a dictionary mapping each state
    dictionary = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6", 7: "7"}

    return dictionary
