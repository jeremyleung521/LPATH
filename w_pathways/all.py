# all.py
#
# Basically a convenient function to run everything!
# Discretize, extract, match, in that order.

from w_pathways import discretize, extract, match
import logging

log = logging.getLogger(__name__)


def main(arguments):
    discretize.main(arguments)
    print(arguments.assign_name)
    extract.main(arguments)
    match.main(arguments)

def entry_point():
    """
    Entry point for discretize, extract, match steps
    """
    import argparse
    from w_pathways import argparser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=argparser.arg_desc)
    argparser.add_discretize_args(parser)
    argparser.add_extract_args(parser)
    argparser.add_match_args(parser)

    args = argparser.process_args(parser)

    discretize.process_assign_args(args)

    log.debug(f'{args}')
    main(args)

if __name__ == "__main__":
    """
    For calling all steps directly. Note all of the parameters are specified manually here.
    """
    import argparse
    args = argparse.Namespace(
        # Discretize Parameters
        input_name="dihedral.npy",  # Input data for state assignment. Something like 'dihedral.npy'.
        output_name="discretized.npy",  # Output file name for the state assignment.
        west_name="multi.h5",  # Name of input HDF5 file (e.g., west.h5)
        assign_name="ANALYSIS/TEST/assign.h5",  # Name of output assign.h5 file
        rcfile="west.cfg",  # west.cfg file
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

        # Extract Parameters
        # Note west_name and assign_name are repeated from above
        source_state_num=0,  # Index of the source state as defined in assign.h5.
        target_state_num=1,  # Index of the target state as defined in assign.h5.
        first_iter=1,  # First iteration to analyze. Inclusive
        last_iter=200,  # Last iteration to analyze. Inclusive. 0 implies it will analyze all labeled iterations.
        trace_basis=True,  # Option to analyze each successful trajectory up till its basis state.
        out_traj=False,  # Option to output trajectory files into `out_dir`. Will take much longer.
        out_traj_ext=".nc",  # Extension of the segment files. Defaults to `seg{out_traj_ext}`.
        out_state_ext="_img.ncrst",  # Extension of the restart files. Defaults to `seg{out_state_ext}`.
        out_top="system.prmtop",  # Name of the parameter file. Name relative to `$WEST_SIM_ROOT/common_files`.
        out_dir="succ_traj",  # Name of directory to output the trajectories.
        hdf5=False,  # Enable if trajectories are saved with the HDF5 Framework in WESTPA.
        rewrite_weights=False,  # Option to zero out the weights of all segments that are not a successful trajectory.
        pcoord=True,  # Option to output the pcoord into the `output.pickle`.
        auxdata=['phi', 'psi'],  # Additional auxiliary data to save into `output.pickle`.
        use_ray=False,  # Enable Ray.
        threads=0,  # How many Ray threads/actors to use. Defaults to 0, which wil use all auto-detected resources.

        # Match Parameters
        # Note west_name, assign_name and out_dir are repeated and removed
        input_pickle='succ_traj/output.pickle',  # Input file name of the pickle from `extract.py`
        dmatrix_remake=True,  # Enable to remake the distance Matrix
        dmatrix_save='distmap.npy',  # If dmatrix_remake is False, load this file instead. Assumed located in {out_dir}.
        dendrogram_threshold=0.5,  # Threshold for the Dendrogram
        dendrogram_show=True,  # Show the Dendrogram using plt.show()
        cl_output='succ_traj/cluster_labels.npy',  # Output path for cluster labels
        file_pattern="west_succ_c{}.h5",  # Pattern to name cluster files
        clusters=None,  # Cluster index to output... otherwise None --> All

        # Others
        debug=False,
    )
    main(args)