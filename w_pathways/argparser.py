# Module that deals with all the argument parsing

import argparse

arg_desc = ''' \
           ARG PARSER
           '''

def check_non_neg(value):
    try:
        value = int(value)
        if value < 0:
            raise argparse.ArgumentTypeError("{} is not a valid iteration number".format(value))
    except ValueError:
        raise Exception("{} is not an integer".format(value))
    return value

def add_args():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description=arg_desc)
    
    iogroup = parser.add_argument_group('input/output options')
    iogroup.add_argument('-W', '--west', '--WEST_H5FILE', '--west-h5file', dest='west_name', default="multi.h5",
        help='''The path to your h5 file. If it's a multi.h5 file from w_multi_west, make sure the --ibstates option successfuly merged your initial and basis states.''')
    iogroup.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name', default='ANALYSIS/TEST/assign.h5', help='')
    iogroup.add_argument('-o', '--output', '--out-traj', '--output-trajectory', dest='out_traj', action='store_true', help='')
    iogroup.add_argument('-oe', '--out-traj-ext', '--output-trajectory-extension', dest='out_traj_ext', default='.nc', type=str, help='')
    iogroup.add_argument('-se', '--out-state-ext', '--output-state-extension', dest='out_state_ext', default='.ncrst', type=str, help='')
    iogroup.add_argument('-ot', '--out-top', '--output-topology', dest='out_top', default='system.prmtop', type=str, help='')
    iogroup.add_argument('-od', '--out-dir', '--output-directory', dest='out_dir', default='succ_traj', type=str, help='')
    iogroup.add_argument('-hdf5', dest='hdf5', action='store_true', help='')
    iogroup.add_argument('-p', '--pcoord', dest='pcoord', action='store_true', help='Output progress coordinate into the pickle file')
    iogroup.add_argument('-a', '--aux', '--AUX', '--auxdata', '--AUXDATA', dest='auxdata', nargs='*', help='''Names of additional auxiliary datasets to be combined''')
    iogroup.add_argument('-aa', '--auxall', action='store_true', help='''Combine all auxiliary datasets. Default: False''')
    iogroup.add_argument('-rw', '--rewrite-weights', action='store_true', help='Copy the H5 files and output individual files where  ')

    parmgroup = parser.add_argument_group('parameters')
    parmgroup.add_argument('-so', '--source', '--source-state', '--SOURCE-STATE', dest='source_state_num', type=check_non_neg, default=0, help='')
    parmgroup.add_argument('-si', '--sink', '--sink-state', '--SINK-STATE', dest='sink_state_num', type=check_non_neg, default=1, help='')
    parmgroup.add_argument('--first', '--first-iter', '--FIRST-ITER', dest='first_iter', type=check_non_neg, default=1, help='')
    parmgroup.add_argument('--last', '--last-iter', '--LAST-ITER', dest='last_iter', type=check_non_neg, default=0, help='')

    opgroup = parser.add_argument_group('Runtime options')
    opgroup.add_argument('-R', '--ray', '--use-ray', dest='use_ray', action='store_true', help='Use Ray work manager.')    
    opgroup.add_argument('-t', '--threads', type=check_non_neg, default=0, help='Number of threads to use with Ray. Default is 0, which uses all available auto-detected resources.')

    args = parser.parse_args()
    
    try:
        import ray
        setattr(args, 'use_ray', True)
    except (ModuleNotFoundError, ImportError):
        pass

    return args
