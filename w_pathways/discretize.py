# discretize.py
#
# Code that allows you to discretize MD trajectories and output as
#
#

import logging
import numpy
import extloader
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
    """
    Main function that executes the whole `match` step. Also called by the
    entry_point() function.

    If it's an HDF5 file, it'll just run w_assign (as ripped from the w_assign west).

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    if arguments.input_name.endswith('.h5'):
        # This basically some logic that's wrapped up in WESTTool.main() for convenience.
        # It needs to be explicitly called like this because the args are captured and set in make_parser_and_process()
        #   which we don't want to call, because we don't want to make a parser.
        #   We just want to process the args that "would've" been captured if called from CLI.
        tool = w_assign.WAssign()

        tool.make_parser_and_process(args=arguments.assign_args)
        # Prepare and instantiate work manager
        tool.wm_env.process_wm_args(arguments.assign_args)
        tool.work_manager = tool.wm_env.make_work_manager()

        tool.process_all_args(args)
        with tool.work_manager:
            if tool.work_manager.is_master:
                tool.go()
            else:
                tool.work_manager.run()
    else:
        input_array = load_file(arguments.input_file)
        out_array = assign(input_array)
        output_file(out_array, arguments.output_file)


def entry_point():
    """
    Entry point for this `match` step.
    """
    import argparse
    from w_pathways import argparser

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=argparser.arg_desc)
    args = argparser.add_discretize_args(parser)

    # Check if any extra w_assign arguments are specified in command line.
    if args.assign_args:
        try:
            import westpa
            from westpa.cli.tools import w_assign
        except ModuleNotFoundError as e:
            print(e)
            raise ModuleNotFoundError("Trying to discretize an HDF5 file but can't import w_assign")

        tool = w_assign.WAssign()
        final_ns = tool.make_parser_and_process(args=args.assign_args)
        setattr(args, 'assign_args', final_ns)
    else:
        default_args = argparse.Namespace(  # These are arguments for w_assign
            verbosity='verbose',  # Verbose or debug
            rcfile='west.cfg',  # west.cfg
            max_queue_length=None,
            we_h5filename='west.h5',  # west.h5 path
            construct_dataset=None,  # If you need some custom auxiliary dataset
            dsspecs=None,
            output='assign.h5',  # Output file
            subsample=None,
            config_from_file=True,  # Read config from rcfile
            scheme='TEST',  # Scheme name
        )
        setattr(args, 'assign_args', default_args)

    log.debug(f'{args}')
    main(args)


if __name__ == "__main__":
    """
    For calling `extract.py` directly. Note all of the parameters are specified manually here.
    """
    import argparse

    args = argparse.Namespace(
        input_name="dihedral.npy",  # Input data for state assignment. Something like 'dihedral.npy'.
        output_name="discretized.npy",  # Output file name for the state assignment.
        west_name="multi.h5",  # Name of input HDF5 file (e.g., west.h5)
        assign_name="ANALYSIS/C7_EQ/assign.h5",  # Name of output assign.h5 file
        rcfile="west.cfg", # west.cfg file
        assign_args=argparse.Namespace(  # These are arguments for w_assign
            verbosity='verbose',  # Verbose or debug
            rcfile="west.cfg",  # west.cfg
            max_queue_length=None,
            we_h5filename='west.h5',  # west.h5 path
            construct_dataset=None,  # If you need some custom auxiliary dataset
            dsspecs=None,
            output='assign.h5',  # Output file
            subsample=None,
            config_from_file=True,  # Read config from rcfile
            scheme='TEST',  # Scheme name
        ),
    )
    main(args)
