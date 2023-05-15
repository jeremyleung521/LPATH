"""
Discretize your MD trajectories (or WE simulations) into states.
"""
import numpy

from mphat.extloader import *
from mphat.io import expanded_load, output_file

log = logging.getLogger(__name__)


def assign(input_array):
    """
    This is an example function for mapping a list of features to state IDs. This should be subclassed
    by passing a similar module (catered to your system) to ``--assign-function``.

    Parameters
    ----------
    input_array : numpy.ndarray or list
        An array generated from expanded_load.

    Returns
    -------
    state_list : list
        A list containing
    """
    state_list = []
    for val in input_array:
        if val[0] >= -180 and val[0] <= -45 and val[1] >= -55 and val[1] <= 30:  # Phi/Psi for Alpha Helix
            state_list.append(0)
        elif val[0] >= 165 and val[0] <= 180 and val[1] >= -55 and val[1] <= 30:
            state_list.append(0)
        elif val[0] >= -170 and val[0] <= -55 and val[1] >= 40 and val[1] <= 100:  # Phi/Psi for C7eq
            state_list.append(1)
        elif val[0] >= 25 and val[0] <= 90 and val[1] >= -55 and val[1] <= 0:  # Phi/Psi for C7ax
            state_list.append(2)
        else:
            state_list.append(3)

    return state_list


def main(arguments):
    """
    Main function that executes the whole ``match`` step.

    If it's been run with ``arguments.we``, it'll just run ``w_assign``.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    if arguments.we:
        try:
            import westpa
            from westpa.cli.tools import w_assign
        except ModuleNotFoundError as e:
            print(e)
            raise ModuleNotFoundError("Trying to discretize an HDF5 file but can't import w_assign")

        if arguments.we and arguments.input_name != arguments.west_name:
            setattr(arguments, 'input_name', arguments.west_name)
            log.debug("Replacing parameter `input_name` with `west_name`")

        if arguments.we and arguments.output_name != arguments.assign_name:
            setattr(arguments, 'output_name', arguments.assign_name)
            log.debug("Replacing parameter `output_name` with `assign_name`")

        if arguments.we and arguments.rcfile != arguments.assign_args.rcfile:
            setattr(arguments, 'rcfile', arguments.assign_args.rcfile)
            log.debug("Replacing parameter `rcfile` with `assign_args.rcfile`")

        # This basically some logic that's wrapped up in WESTTool.main() for convenience.
        # It needs to be explicitly called like this because the args are captured and set in make_parser_and_process()
        #   which we don't want to call, because we don't want to make a parser.
        #   We just want to process the args that "would've" been captured if called from CLI.

        tool = w_assign.WAssign()

        # Prepare and instantiate work manager
        tool.wm_env.process_wm_args(arguments.assign_args)
        tool.work_manager = tool.wm_env.make_work_manager()

        tool.process_all_args(arguments.assign_args)
        with tool.work_manager:
            if tool.work_manager.is_master:
                tool.go()
            else:
                tool.work_manager.run()
    else:
        input_array = expanded_load(arguments.input_name, arguments.stride)

        # Replacing assign_func with what's given
        if arguments.assign_func != 'default_assign':
            import sys
            import os
            sys.path.append(os.getcwd())

            assign = get_object(arguments.assign_func)
            log.info(f'INFO: Replaced assign() with {arguments.assign_func}')

        out_array = assign(input_array)

        # Warning
        n_states = len(set(out_array))
        if n_states < 3:
            log.info(f'Only {n_states} defined, including the "unknown" state. This should be fine for \
                       ``mphat extract`` purposes but will likely produce bad quality pattern matching results \
                       further downstream. Please consider running ``mphat discretize`` with more states or \
                       plan to reassign the trajectory using the ``reassign-method`` option in ``mphat match``.')

        # Output
        output_file(out_array, arguments.output_name)
