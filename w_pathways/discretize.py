# discretize.py
#
# Code that allows you to discretize MD trajectories and output as
#
#

import logging
import numpy
import w_pathways
from tqdm.auto import tqdm

log = logging.getLogger(__name__)


def assign(input_array):
    state_list = []
    for val in tqdm(input_array):
        if val[0] >= -180 and val[0] <= -45 and val[1] >= -55 and val[1] <= 30:  # Phi/Psi for Alpha Helix
            state_list.append(0)
        elif val[0] >= 165 and val[0] <= 180 and val[1] >= -55 and val[1] <= 30:
            state_list.append(0)
        elif val[0] >= -170 and val[0] <= -55 and val[1] >= 40 and val[1] <= 100:  # Phi/Psi for C7eq
            state_list.append(1)
        elif val[0] >= 25 and val[0] <= 90 and val[1] >= -55 and val[1] <= 0:  # Phi/Psi for C7ax
            state_list.append(2)
        else:
            state_list.append(-1)

    return state_list


def load_file(input_file):
    """
    Parameters
    ----------
    input_file: str
        Path of the data file to be used to assign states.

    Returns
    -------
    data: numpy.array
        A numpy array of data used to assign states.
    """
    if input_file.endswith('.npy'):
        data = numpy.load(input_file)
    else:
        data = numpy.loadtxt(input_file)
        # data = numpy.loadtxt(input_file, usecols=(1,2), skiprows=1)
    return data


def output_file(out_array, output_name):
    n = numpy.asarray(out_array)
    numpy.save(output_name, n)


def main(arguments):
    input_array = load_file(arguments.input_file)
    out_array = assign(input_array)
    output_file(out_array, arguments.output_file)


def entry_point():
    import argparse
    from w_pathways import argparser

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=argparser.arg_desc)
    args = argparser.add_discretize_args(parser)
    log.debug(f'{args}')
    main(args)


if __name__ == "__main__":
    # If calling discretize.py directly...
    import argparse
    args = argparse.Namespace(
        input_name="dihedral.npy",
        output_name="discretized.npy",
    )
    main(args)
