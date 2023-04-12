'''
This module deals with all the argument parsing.
'''
import argparse
import logging

log = logging.getLogger(__name__)

arg_desc = '''
           mPHAT: minimal Pathway Analysis Histogram Analysis of Trajectories
           ==================================================================
           '''


def check_non_neg(value):
    try:
        value = int(value)
        if value < 0:
            raise argparse.ArgumentTypeError("{} is not a valid iteration number".format(value))
    except ValueError:
        raise Exception("{} must be an integer.".format(value))
    return value


def add_discretize_args(parser):
    """
    This block process all the necessary arguments for the `discretize.py` module.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        A parser passed in from each tool. Separated from each function because the
        catch-all tool to run everything in succession will only have 1 parser.
    """
    discrete_io = parser.add_argument_group('Discretize input/output options')
    discrete_io.add_argument('-I', '--input', dest='input_name', default="west.h5",
                             help='''The path to your input file for discretization. If it's a west.h5 file from 
                             WESTPA, this will automatically run w_assign. Else, it would assume it's a text file of 
                             features.''')
    discrete_io.add_argument('-O', '--output', dest='output_name', default="states.npy",
                             help='''The path to your output numpy file for after discretization. ''')

    discrete_io.add_argument('-ar', '--assign-args', '--assign-arguments', dest='assign_args', type=str, default='',
                             help='''A string of arguments to pass onto w_assign as you would imput in the command 
                             line to `w_assign`. Either use the defaults (leave blank) or at a minimum, you need to add 
                             in `--states-from-config --scheme NAME_OF_SCHEME` to read config from your `west.cfg`''')

    try:
        discrete_io.add_argument('-W', '--west', '--WEST_H5FILE', '--west-h5file', dest='west_name', default="west.h5",
                                 help='''The path to your h5 file. If it's a multi.h5 file from w_multi_west, make sure 
                                 the `--ibstates` option successfully merged your initial and basis states.''')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        discrete_io.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name',
                                 default='ANALYSIS/TEST/assign.h5', help='')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        discrete_io.add_argument('-r', '--rcfile', metavar='RCFILE', dest='rcfile', default='west.cfg',
                                 help='use RCFILE as the WEST run-time configuration file (default: %(default)s)')
    except argparse.ArgumentError as e:
        log.debug(e)

    try:
        discrete_io.add_argument('--debug', action='store_true', help='Enable debug mode.')
    except argparse.ArgumentError as e:
        log.debug(e)


def add_extract_args(parser):
    """
    This block process all the necessary arguments for the "extract.py" module.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        A parser passed in from each tool. Separated from each function because the
        catch-all tool to run everything in succession will only have 1 parser.
    """
    iogroup = parser.add_argument_group('Extract input/output options')
    try:
        iogroup.add_argument('-W', '--west', '--WEST_H5FILE', '--west-h5file', dest='west_name', default="west.h5",
                             help='''The path to your h5 file. If it's a multi.h5 file from w_multi_west, make sure 
                             the `--ibstates option successfully merged your initial and basis states.''')
        iogroup.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name',
                             default='ANALYSIS/TEST/assign.h5', help='')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        iogroup.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name',
                             default='ANALYSIS/TEST/assign.h5', help='')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        iogroup.add_argument('-od', '--out-dir', '--output-directory', dest='out_dir', default='succ_traj',
                             type=str, help='')
    except argparse.ArgumentError as e:
        log.debug(e)

    iogroup.add_argument('-oj', '--out-traj', '--output-trajectory', dest='out_traj',
                         action='store_true', help='')
    iogroup.add_argument('-oe', '--out-traj-ext', '--output-trajectory-extension', dest='out_traj_ext',
                         default='.nc', type=str, help='')
    iogroup.add_argument('-se', '--out-state-ext', '--output-state-extension', dest='out_state_ext',
                         default='.ncrst', type=str, help='')
    iogroup.add_argument('-ot', '--out-top', '--output-topology', dest='out_top', default='system.prmtop',
                         type=str, help='')

    iogroup.add_argument('-hdf5', '--hdf5', dest='hdf5', action='store_true', help='')

    parmgroup = parser.add_argument_group('Extract parameters')
    parmgroup.add_argument('-ss', '--source', '--source-state', '--SOURCE-STATE', dest='source_state_num',
                           type=check_non_neg, default=0, help='')
    parmgroup.add_argument('-ts', '--target', '--target-state', '--SINK-STATE', dest='target_state_num',
                           type=check_non_neg, default=1, help='')
    parmgroup.add_argument('--first', '--first-iter', '--FIRST-ITER', dest='first_iter', type=check_non_neg,
                           default=1, help='')
    parmgroup.add_argument('--last', '--last-iter', '--LAST-ITER', dest='last_iter', type=check_non_neg,
                           default=0, help='')
    parmgroup.add_argument('--trace-basis', '-b', dest='trace_basis', action='store_true', help='')
    parmgroup.add_argument('--pcoord', '-p', dest='pcoord', action='store_true',
                           help='Output progress coordinate into the pickle file')
    parmgroup.add_argument('-a', '--aux', '--AUX', '--auxdata', '--AUXDATA', dest='auxdata', nargs='*',
                           help='''Names of additional auxiliary datasets to be combined''')
    parmgroup.add_argument('-aa', '--auxall', action='store_true',
                           help='''Combine all auxiliary datasets. Default: False''')
    parmgroup.add_argument('--rewrite-weights', '-rw', action='store_true',
                           help='Copy the H5 files and output individual files where  ')

    opgroup = parser.add_argument_group('Extract runtime options')
    opgroup.add_argument('--use-ray', '-R', '--ray', dest='use_ray', action='store_true', help='Use Ray work manager.')
    opgroup.add_argument('--no-ray', '-NR', dest='no_ray', action='store_true',
                         help='Do not use Ray. This overrides --use-ray.')
    opgroup.add_argument('-t', '--threads', type=check_non_neg, default=0, help='Number of threads to use \
     with Ray. Default is 0, which uses all available auto-detected resources.')
    try:
        opgroup.add_argument('--debug', action='store_true', help='Enable debug mode.')
    except argparse.ArgumentError as e:
        log.debug(e)


def add_match_args(parser):
    """
    This block process all the necessary arguments for the "match.py" module.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        A parser passed in from each tool. Separated from each function because the
        catch-all tool to run everything in succession will only have 1 parser.
    """
    match_io = parser.add_argument_group('Match input/output options')
    try:
        match_io.add_argument('-W', '--west', '--WEST_H5FILE', '--west-h5file', dest='west_name', default='west.h5',
                              help='''The path to your h5 file. If it's a multi.h5 file from w_multi_west, make sure 
                              the `--ibstates option successfully merged your initial and basis states.''')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        match_io.add_argument('-A', '--assign', '--assign-h5file', '--ASSIGN-H5FILE', dest='assign_name',
                              default='ANALYSIS/TEST/assign.h5', help='')
    except argparse.ArgumentError as e:
        log.debug(e)
    try:
        match_io.add_argument('-od', '--out-dir', '--output-directory', dest='out_dir', default='succ_traj',
                              type=str, help='')
    except argparse.ArgumentError as e:
        log.debug(e)
    match_io.add_argument('--pickle', '--input-pickle', dest='input_pickle', default='succ_traj/output.pickle',
                          type=str, help='Path to pickle object from `extract`')
    match_io.add_argument('-cd', '--cl-out-dir', '--cluster-label-output-directory', dest='cl_output',
                          default='succ_traj', type=str, help='')
    match_io.add_argument('-fp', '--fp', '--file-pattern', dest='file_pattern',
                          default="west_succ_c{}.h5", type=str, help='Pattern to name cluster files.')
    match_io.add_argument('-ex', '--fp', '--export-h5', dest='export_h5',
                          action='store_true', help='Export each cluster as an independent H5 file.')

    match_parm = parser.add_argument_group('Match parameters')
    match_parm.add_argument('-dt', '--dendro-threshold', '--dendrogram-threshold', dest='dendrogram_threshold',
                            type=check_non_neg, default=0.5,
                            help='Horizontal threshold line for dendrogram.')
    match_parm.add_argument('-ds', '--dendro-show', '--dendrogram-show', dest='dendrogram_show', action='store_true',
                            help='Show dendrogram with `plt.show()`.')
    match_parm.add_argument('-dh', '--dendro-hide', '--dendrogram-hide', dest='dendrogram_hide', action='store_true',
                            help='Do not show dendrogram. Overrides `--dendrogram-show`.')
    match_parm.add_argument('-c', '--clusters', dest='clusters', default=None, nargs='*',
                            help='Clusters to export. 0-indexed. Default: None')

    match_op = parser.add_argument_group('Match runtime options')
    match_op.add_argument('--reassign', '-ra', dest='reassign_method', default='reassign_identity', type=str,
                          help='Reassign Method to use.')
    match_op.add_argument('--remake', '-dR', dest='dmatrix_remake', action='store_false',
                          help='Remake distance matrix.')
    match_op.add_argument('--remake-file', '-dF', dest='dmatrix_save', type=str, default='distmap.npy',
                          help='Path to pre-calculated distance matrix. Assumed to be in `out_dir`.')
    match_op.add_argument('--no-remake', '-nR', dest='no_remake', action='store_true',
                          help='Do not remake distance matrix. This overrides `--remake.`')

    try:
        match_op.add_argument('--debug', action='store_true', help='Enable debug mode.')
    except argparse.ArgumentError as e:
        log.debug(e)


def process_args(parser):
    args = parser.parse_args()

    # Automatically turn on Ray unless no_ray is specified.
    try:
        import ray
        if args.no_ray is False:
            setattr(args, 'use_ray', True)
        elif args.no_ray is True and args.use_ray is True:
            setattr(args, 'use_ray', False)
    except (ModuleNotFoundError, ImportError, AttributeError):
        pass

    # Turn Debugging on!
    if args.debug is True:
        logging.basicConfig(level=logging.DEBUG)

    return args
