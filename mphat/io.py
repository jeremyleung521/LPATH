"""
I/O Operations
"""

import numpy

def load_file(input_file, stride):
    """
    Parameters
    ----------
    input_file: str
        Path of the data file to be used to assign states.

    stride : int
        Dictates how often to load in the data. Only used in Standard MD.

    Returns
    -------
    data: numpy.array
        A numpy array of data used to assign states.

    """
    if input_file.endswith('.npy'):
        data = numpy.load(input_file)[::stride]
    else:
        data = numpy.loadtxt(input_file)[::stride]
        # data = numpy.loadtxt(input_file, usecols=(1,2), skiprows=1)

    return data


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