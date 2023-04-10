# Module that deals with all the argument parsing

import argparse
import logging

arg_desc = ''' \
           ARG PARSER
           '''

def check_non_neg(value):
    try:
        value = int(value)
        if value < 0:
            raise argparse.ArgumentTypeError("{} is not a valid iteration number".format(value))
    except ValueError:
        raise Exception("{} must be an integer.".format(value))
    return value

class RC:
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=arg_desc)

    def add_discretize_args(self):
        """
        This block process all the necessary for the "discretize.py" module.

        Returns
        -------
        args: NameSpace
            Object with all of the necessary arguments.

        """
        discrete_io = self.parser.add_argument_group('discretize input/output options')
        discrete_io.add_argument('-I', '--input', dest='input_name', default="west.h5",
            help='''The path to your input file for discretization. If it's a west.h5 file from WESTPA, 
            this will automatically run w_assign. Else, it would assume it's a text file of "features"''')
        discrete_io.add_argument('-O', '--output', dest='output_name', default="states.npy",
            help='''The path to your output numpy file for after discretization. ''')

        return args

    def add_extract_args(self):
        """
        This block process all the necessary for the "extract.py" module.

        Returns
        -------
        args: NameSpace
            Object with all of the necessary arguments.

        """
        iogroup = self.parser.add_argument_group('input/output options')
        iogroup.add_argument('-W', '--west', '--WEST_H5FILE', '--west-h5file', dest='west_name', default="multi.h5",
            help='''The path to your h5 file. If it's a multi.h5 file from w_multi_west, make sure the \
            --ibstates option successfully merged your initial and basis states.''')
        iogroup.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name',
                             default='ANALYSIS/TEST/assign.h5', help='')
        iogroup.add_argument('--output', '-o', '--out-traj', '--output-trajectory', dest='out_traj',
                             action='store_true', help='')
        iogroup.add_argument('-oe', '--out-traj-ext', '--output-trajectory-extension', dest='out_traj_ext',
                             default='.nc', type=str, help='')
        iogroup.add_argument('-se', '--out-state-ext', '--output-state-extension', dest='out_state_ext',
                             default='.ncrst', type=str, help='')
        iogroup.add_argument('-ot', '--out-top', '--output-topology', dest='out_top', default='system.prmtop',
                             type=str, help='')
        iogroup.add_argument('-od', '--out-dir', '--output-directory', dest='out_dir', default='succ_traj',
                             type=str, help='')
        iogroup.add_argument('-hdf5', '--hdf5', dest='hdf5', action='store_true', help='')
        iogroup.add_argument('--pcoord', '-p', dest='pcoord', action='store_true',
                             help='Output progress coordinate into the pickle file')
        iogroup.add_argument('-a', '--aux', '--AUX', '--auxdata', '--AUXDATA', dest='auxdata', nargs='*',
                             help='''Names of additional auxiliary datasets to be combined''')
        iogroup.add_argument('-aa', '--auxall', action='store_true',
                             help='''Combine all auxiliary datasets. Default: False''')
        iogroup.add_argument('--rewrite-weights', '-rw', action='store_true',
                             help='Copy the H5 files and output individual files where  ')

        parmgroup = self.parser.add_argument_group('parameters')
        parmgroup.add_argument('-ss', '--source', '--source-state', '--SOURCE-STATE', dest='source_state_num',
                               type=check_non_neg, default=0, help='')
        parmgroup.add_argument('-ts', '--target', '--target-state', '--SINK-STATE', dest='target_state_num',
                               type=check_non_neg, default=1, help='')
        parmgroup.add_argument('--first', '--first-iter', '--FIRST-ITER', dest='first_iter', type=check_non_neg,
                               default=1, help='')
        parmgroup.add_argument('--last', '--last-iter', '--LAST-ITER', dest='last_iter', type=check_non_neg,
                               default=0, help='')
        parmgroup.add_argument('--trace-basis', '-b', dest='trace_basis', action='store_true', help='')

        opgroup = self.parser.add_argument_group('Runtime options')
        opgroup.add_argument('--use-ray', '-R', '--ray', dest='use_ray', action='store_true', help='Use Ray work manager.')
        opgroup.add_argument('--no-ray', '-NR', dest='no_ray', action='store_true',
                             help='Do not use Ray. This overrides --use-ray')
        opgroup.add_argument('-t', '--threads', type=check_non_neg, default=0, help='Number of threads to use \
         with Ray. Default is 0, which uses all available auto-detected resources.')
        opgroup.add_argument('--debug', action='store_true', help='Enable debug mode.')


        args = self.parser.parse_args()

        # Automatically turn on Ray unless
        try:
            import ray
            if args.no_ray is False:
                setattr(args, 'use_ray', True)
            elif args.no_ray is True and args.use_ray is True:
                setattr(args, 'use_ray', False)
        except (ModuleNotFoundError, ImportError):
            pass

        # Turn Debugging on!
        if args.debug is True:
            logging.basicConfig(level=logging.DEBUG)

        return args

