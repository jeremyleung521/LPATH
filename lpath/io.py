"""
I/O Operations
"""
import numpy
from ast import literal_eval
from os import mkdir
from ._logger import Logger

log = Logger().get_logger(__name__)

default_dendrogram_colors = ['tomato', 'dodgerblue', 'orchid', 'mediumseagreen', 'darkorange', 'mediumpurple', 'grey']


def load_file(input_file, stride=1):
    """
    The actual loading function.

    Parameters
    ----------
    input_file: str
        Path (or lists of paths) of the data file used to assign states.

    stride: int, default : 1
        Step to stride file. Defaults to 1, which returns everything.

    Returns
    -------
    data: numpy.array
        A numpy array of data used to assign states.

    """
    try:
        # Treat it as text first.
        data = numpy.loadtxt(input_file)[::stride]
        # data = numpy.loadtxt(input_file, usecols=(1,2), skiprows=1)
    except UnicodeDecodeError:
        # Assumes it's a binary file otherwise...
        log.debug(f'{input_file} is not a unicode text file. Attempting to load file like a binary file.')
        data = numpy.load(input_file, allow_pickle=True)[::stride]

    return data


def expanded_load(input_str, stride):
    """
    The expanded loading function that actually deals with lists. Attempts to literal_eval
    it before processing.

    Parameters
    ----------
    input_str : string or list of strings
        Path or a list of paths of files to be opened. Files will be step-sliced
        (using the ``stride`` parameter) independently.

    stride : int
        Dictates how often to load in the data. Only used in Standard MD.

    Returns
    -------
    final_object : list or numpy.ndarray
        The loaded items. If input was a list of files, they would be sliced independently
        and then concatenated together.

    """
    try:
        expanded_files = literal_eval(input_str)
        if not isinstance(expanded_files, list):
            raise ValueError
    except ValueError:
        final_object = load_file(input_str, stride)
        return final_object

    # Loop through all strings in expanded_files since it's a list for sure.
    try:
        final_object = load_file(expanded_files[0], stride)
        for file in expanded_files[1:]:
            final_object = numpy.append(final_object, load_file(file, stride))
        return final_object
    except FileNotFoundError as e:
        raise FileNotFoundError(f'Could not open files for processing: {e.args}')


def output_file(out_array, output_name):
    """
    Function to output an array.

    Parameters
    ----------
    out_array: numpy.ndarray
        Array to be outputted.

    output_name: str
        Name of the output file.

    """
    n = numpy.asarray(out_array)
    numpy.save(output_name, n)


def make_dir(args):
    """
    Make directories according to the argparse.Namespace object.

    Parameters
    ----------
    args : argparse.Namespace
        Passed directly from argparser.

    """
    if args.step_name in ['extract', 'all']:
        try:
            mkdir(args.out_dir)
        except FileExistsError:
            log.warning(f"Folder ``{args.out_dir}`` already exists. Files within might be overwritten.")

    if args.step_name in ['match', 'plot', 'all']:
        try:
            mkdir(args.out_path)
        except FileExistsError:
            log.warning(f"Folder ``{args.out_path}`` already exists. Files within might be overwritten.")


class EmptyOutputError(Exception):
    """
    Custom Error for cases when extract has empty output.

    """
    def __init__(self, message="No successful trajectories extracted."):
        self.message = message
        super().__init__(self.message)
