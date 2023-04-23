"""
Extract successful trajectories from MD trajectories (or WE simulations).
"""
# Code that traces through an assign.h5 to generate a list of lists containing
#     successful source -> target transitions. 
#
# This version has ray parallelization.
#
# Modified to work with pathway analysis code. 
#
# Main function includes:
#   * Output a list containing traces of source -> sink iter/seg numbers
#   * Output the trajectories of the traces
#   * Rewrite a H5 file so any non-successful trajectories have 0 weight.
#
# Defaults to the first passage (from basis state to first touch target)
#     for all of the above.
# Make `trace_basis` False to trace only the barrier crossing time.

import logging
import pickle
from collections import Counter
from copy import deepcopy
from os import mkdir
from shutil import copyfile

import numpy
from tqdm.auto import trange

from mphat.io import load_file, expanded_load

log = logging.getLogger(__name__)


# Here are functions for standard MD.
def find_min_distance(ref_value, indices):
    """
    Search for the closest value in indices that is >= to ref_value.

    Parameters
    ----------
    ref_value : int or float
        Reference point you want to anchor the search.

    indices : list or numpy.ndarray
        A list or array of potential values you want to search for.

    Returns
    -------
    minimum value : float, int
        The closest value to ref_value in indices.

    """
    # TODO: SPEED THIS UP
    # THIS ACTUALLY TAKES QUITE SOME TIME.

    return min([v for v in indices if v >= ref_value] or [None])


def clean_self_to_self(input_array):
    """
    Clean up duplicates which might contain self to self transitions.

    Parameters
    ----------
    input_array : list
        A list or numpy array of the shape (n_transitions, 2).

    Returns
    -------
    output_array : numpy.ndarray
        A reduced list or numpy array of the shape (n_transitions, 2).

    """
    output_array = numpy.asarray(input_array)
    full_count = Counter(output_array[:, 1])

    # Excluding those that appear only once
    reduced_keys = [key for (key, count) in full_count.items() if count > 1]

    # For debugging purposes.
    # log.debug(f"Full Count: {full_count}")
    # log.debug(f"Reduced Keys: {reduced_keys}")
    # Running backwards so indices are maintained.
    for delete in reduced_keys[::-1]:
        # Determine indices of where the duplicates happen
        pop_list = numpy.argwhere(output_array[:, 1] == delete).flatten()
        pop_list.sort()
        # log.debug(f"All indices where {delete} occur: {pop_list}")
        # Remove them except for the last instance
        for j in pop_list[::-1][1:]:
            input_array.pop(j)

    output_array = numpy.asarray(input_array)

    return output_array


def count_tmatrix_row(source_index, trajectory, n_states, target_num):
    """
    Count transitions for the source --> states row for the weights. Used to
    calculate the weights of each successful trajectory.

    Parameters
    ----------
    source_index : numpy.ndarray
        An array of indices where it last visited the source state.

    trajectory : numpy.ndarray
        A list of the states as inputted.

    n_states : int
        Number of total states defined. Does not include the Unknown State.

    source_num : int
        Index of the source state as defined in ``discretize``.

    target_num : int
        Index of the target state as defined in ``discretize``.

    Returns
    -------
    st_weight : float
        Total weight of all the source --> target transitions.

    """
    count_row = numpy.zeros(n_states)
    for istate in source_index:
        for jstate in trajectory[istate + 1:]:
            # If it isn't in the Unknown State.
            if jstate != n_states:
                count_row[jstate] += 1
                break

    # Row Normalize for the probability
    count_row /= sum(count_row)
    log.debug(f'Source row for T-matrix: {count_row}')
    st_weight = count_row[target_num]

    return st_weight


def find_transitions(input_array, source_index, target_index):
    """
    Find all successful transitions in standard MD simulations.

    Parameters
    ----------
    input_array : numpy.ndarray
        An array of states assignments. Should be of shape (n_frames).

    source_index : int or float
        The assignment of the source state.

    target_index : int or float
        The assignment of the target state.

    Returns
    -------
    source_indices : numpy.ndarray
        An array of all indices where the input data visited the source state.

    target_indices : numpy.ndarray
        An array of all indices where the input data visited the target state.

    transitions : numpy.ndarray
        An array of shape (n_transitions, 2) showing all steps.

    """
    source_indices = numpy.argwhere(input_array == source_index).flatten()
    target_indices = numpy.argwhere(input_array == target_index).flatten()

    log.debug(f'Indices of source state {source_index}: {source_indices}')
    log.debug(f'Indices of target state {target_index}: {target_indices}')

    # Doing some sanity checks...
    assert len(source_indices) > 0, f"None of your frames were assigned to the source state {source_index}."
    assert len(target_indices) > 0, f"None of your frames were assigned to the target state {target_index}."
    # This following check is taking too long.
    # assert any(target > source for target in target_indices for source in source_indices), \
    #     f"None of your target state assignments {target_index} occur after the \
    #        last visit to the source state {source_index}."

    # Now do the calculations
    transitions = []
    for val in source_indices:
        check = find_min_distance(val, target_indices)
        if check is not None:
            transitions.append([val, check])

    log.debug(f'Indices of all successful transitions: {transitions}')

    return source_indices, target_indices, transitions


def create_pickle_obj(transitions, states, weight, features=None):
    """
    Main function that transforms a list of frame transitions into the pickle object. For standard simulations only.

    Parameters
    ----------
    transitions : list
        A list of shape (successful transitions, 2). Indicates the start/end frame a transition
        has been made.

    states : list
        A list of the states as inputted.

    weight : float
        Weight of each successful transition.

    features : list or numpy.ndarray or None
        If specified by the user, you can save extra information in to the pickle object.

    Returns
    -------
    output_list : list
        An list to be outputted, prepared for the ``output.pickle``.

    """
    if features is None:
        ad_arr = []
    else:
        if isinstance(features, numpy.ndarray):
            ad_arr = features.tolist()
        else:
            ad_arr = list(features)

    output_list = []

    # Going through each transition
    for idx, transition in enumerate(transitions):
        indv_trace = []
        # Add all the frames between target + source. This is controlled by stride, which dictates how often you load.
        for j in range(int(transition[0]), int(transition[1] + 1)):
            indv_trace.append([1, idx, states[j], *ad_arr[j], weight])
        output_list.append(deepcopy(indv_trace))

    return output_list


def standard(arguments):
    """
    Main function that executes the whole standard `extract` procedure.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    input_array = load_file(arguments.extract_input, arguments.stride)
    n_states = len(input_array)

    if arguments.pcoord is True:
        if arguments.featurization_name is None:
            log.warning('WARNING: The --pcoord flag is specified but no file is \
                         specified with --extract-pcoord. Skipping output.')
            features = None
        else:
            features = expanded_load(arguments.featurization_name, arguments.stride)
    else:
        features = None

    source_index, target_index, transitions = find_transitions(input_array, arguments.source_state_num,
                                                               arguments.target_state_num)

    new_transitions = clean_self_to_self(transitions)

    weight = count_tmatrix_row(source_index, input_array, n_states, arguments.target_state_num)

    # Generate and write pickle object.
    final_obj = create_pickle_obj(transitions, input_array, weight / len(transitions), features)

    with open(f"{arguments.out_dir}/{arguments.extract_output}", "wb") as fo:
        pickle.dump(final_obj, fo)

    with open(f"{arguments.out_dir}/frame_info.pickle", "wb") as fo:
        pickle.dump(new_transitions, fo)


# Here are functions for WE.
def we(arguments):
    """
    Main function that executes the whole WE `extract` procedure.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    import h5py
    import westpa.analysis as wa

    def process_ad_arr(pcoord, auxdata, stride, iwalker):
        """
        Inner function that deals with preparing ad_arr

        Parameters
        ----------
        pcoord : bool
            Marks whether to add progress coordinate in ad_arr or not.

        auxdata : list or None
            A list of auxiliary datasets to add or None

        stride : int
            How much to stride the current data.


        iwalker : westpa.analysis.core.Walker
            A walker object from ``westpa.analysis`` of relevant segment.

        Returns
        -------
        ad_arr : list
            A list of relevant data to be outputted.

        """
        ad_arr = []
        total_frames = iwalker.pcoords.shape[0]
        stride_step = total_frames // stride
        if pcoord is True:
            if len(iwalker.pcoords.shape) > 1:
                # for item in iwalker.pcoords[::-stride_step]:
                for item in iwalker.pcoords[-1]:
                    ad_arr.append(item)
            else:
                # ad_arr.append(iwalker.pcoords[::-stride_step])
                ad_arr.append(iwalker.pcoords[-1])

        if auxdata is not None:
            if len(auxdata) == 0:
                # auxdata is set to 0 when auxall is called, so grabbing all the aux datasets.
                auxdata = list(iwalker.auxiliary_data.keys())
            for dataset_name in auxdata:
                if len(iwalker.auxiliary_data[dataset_name].shape) > 1:
                    # for item in iwalker.auxiliary_data[dataset_name][::-stride_step]:
                    for item in iwalker.auxiliary_data[dataset_name][-1]:
                        ad_arr.append(item)
                else:
                    # ad_arr.append(iwalker.auxiliary_data[dataset_name][::-stride_step])
                    ad_arr.append(iwalker.auxiliary_data[dataset_name][-1])

        return ad_arr, range(total_frames, -1, stride_step)

    # Defining functions if we're using ray...
    if arguments.use_ray:
        import ray

        @ray.remote
        class TraceActor:
            """
            An actor class that has its own copy of west.h5 file and assign.h5 file.

            """
            def __init__(self, assign_name, run_name):
                self.assign_file = h5py.File(assign_name, 'r')
                self.h5run = wa.Run(run_name)

            def trace_seg_to_last_state(
                    self,
                    source_state_num,
                    target_state_num,
                    iteration_num,
                    segment_num,
                    first_iter,
                    trace_basis,
                    out_traj,
                    out_traj_ext,
                    out_state_ext,
                    out_top,
                    hdf5,
                    pcoord,
                    auxdata,
                    stride,
            ):
                """
                Code that traces a seg to frame it leaves source state. Can run to export trajectories too!

                Returns
                =======

                indv_trace : lst of lst
                    A list of lists containing the iteration and segment numbers of the trace. None if the
                    trajectory does not end in the target state.

                indv_traj : obj
                    A BasicMDTrajectory() or HDF5MDTrajectory() object or None. Basically a subclass of the
                    MDTraj Trajectory object.

                frame_info : lst of lst
                    A list of list containing the frame number of each iteration. Goes backwards from the last frame.
                    Returns None if it does not end in the target state.

                """
                run = self.h5run
                assign_file = self.assign_file
                indv_trace = []
                trace = run.iteration(iteration_num).walker(segment_num).trace()
                traj_len = (
                        len(run.iteration(iteration_num).walker(segment_num).pcoords) - 1
                )  # Length is number of frames in traj + 1 (parent); only caring about the number of frames

                # Going through segs in reverse order
                for iwalker in reversed(trace):
                    ad_arr, stride_steps = process_ad_arr(pcoord, auxdata, stride, iwalker)
                    ad_arr = numpy.asarray(ad_arr)

                    weight = iwalker.weight
                    corr_assign = assign_file["statelabels"][
                        iwalker.iteration.summary.name - 1, iwalker.segment_summary.name
                    ]
                    if iwalker.iteration.summary.name == iteration_num:
                        # Dealing with cases where this is the first iteration we're looking at
                        # If multiple, taking only the first instance we reached tstate
                        term_frame_num = numpy.where(corr_assign == target_state_num)[0][0]
                        if trace_basis is False:
                            if source_state_num in corr_assign[:term_frame_num + 1]:
                                # Went from source to target in one iteration. neat.
                                source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                                indv_trace.append([iteration_num, segment_num, corr_assign[-1], *ad_arr, weight])
                                break
                        # Just a normal iteration where we reached target state and that's it.
                        indv_trace.append([iteration_num, segment_num, corr_assign[-1], *ad_arr, weight])

                    elif iwalker.iteration.summary.name != iteration_num:
                        # Dealing with cases where we're in other iterations
                        if source_state_num not in corr_assign:
                            # The segment did not visit the source this iteration.
                            if iwalker.iteration.summary.name != iteration_num:
                                if target_state_num in corr_assign:
                                    # Hey, this is a target > target transition.
                                    # Breaking out...
                                    return None, None, None
                                else:
                                    # If traj hasn't been in state A, and not in target iteration...
                                    # add the whole iteration into list.
                                    # Also making sure it's not a target -> target transition
                                    indv_trace.append(
                                        [iwalker.iteration.summary.name, iwalker.segment_summary.name, corr_assign[-1],
                                         *ad_arr, weight]
                                    )
                        else:
                            # This else captures cases where it was in the source in this iteration
                            # Looking for the last frame it was in source state
                            source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                            if target_state_num in corr_assign[source_frame_num:]:
                                # Catching the ouchie case where it went from source -> target -> target. Woopsies!
                                return None, None, None
                            else:
                                # Final case where it's definitely source -> target
                                indv_trace.append(
                                    [iwalker.iteration.summary.name, iwalker.segment_summary.name, corr_assign[-1],
                                     *ad_arr, weight]
                                )
                                if trace_basis is False:
                                    break

                frame_info = []
                try:
                    source_frame_num
                except NameError:
                    source_frame_num = 0

                # Total number of frames necessary
                start_trace = (traj_len - source_frame_num) + ((len(indv_trace) - 1) * traj_len)
                end_trace = traj_len - term_frame_num
                frame_info.append(term_frame_num)
                for i in range(len(indv_trace) - 1, 1, -1):
                    frame_info.append(traj_len)
                frame_info.append(source_frame_num)

                # Block for outputting the traj
                if out_traj:
                    if hdf5 is False:
                        trajectory = wa.BasicMDTrajectory(
                            traj_ext=out_traj_ext, state_ext=out_state_ext, top=out_top
                        )
                    elif hdf5 is True:
                        trajectory = wa.HDF5MDTrajectory()
                    else:
                        log.warning("unable to output trajectory")
                        return indv_trace, None, frame_info

                    # Extra check such that it won't output past first_iter.
                    if first_iter != 1:
                        max_length = iteration_num - first_iter + 1
                        trace = run.iteration(iteration_num).walker(segment_num).trace(max_length=max_length)

                    indv_traj = trajectory(trace)
                    if trace_basis is False:
                        indv_traj = indv_traj[-start_trace:-end_trace]
                else:
                    indv_traj = None

                # print(indv_trace)

                return indv_trace, indv_traj, frame_info
    else:
        # Or the case where we won't!
        def trace_seg_to_last_state(
                source_state_num,
                target_state_num,
                new_file,
                assign_file,
                iteration_num,
                segment_num,
                trace_basis,
                out_traj,
                out_traj_ext,
                out_state_ext,
                out_top,
                hdf5,
                tqdm_bar,
                pcoord,
                auxdata,
                stride,
        ):
            """
            Code that traces a seg to frame it leaves source state. Can run to export trajectories too!
            Returns
            =======
            indv_trace : lst of lst
                A list of list containing the iteration and segment numbers of the trace. None if the trajectory
                does not end in the target state.
            indv_traj : obj
                A BasicMDTrajectory() or HDF5MDTrajectory() object or None. Basically a subclass of the
                MDTraj Trajectory object.
            frame_info : lst of lst
                A list of list containing the frame number of each iteration. Goes backwards from the last frame.
                Returns None if traj does not end in the target state.
            """
            run = wa.Run(new_file)
            indv_trace = []
            trace = run.iteration(iteration_num).walker(segment_num).trace()
            traj_len = (
                    len(run.iteration(iteration_num).walker(segment_num).pcoords) - 1
            )  # Length is number of frames in traj + 1 (parent); only caring about the number of frames

            tqdm_bar.set_description(f"tracing {iteration_num}.{segment_num}")
            # Going through segs in reverse order
            for iwalker in reversed(trace):
                ad_arr, stride_steps = process_ad_arr(pcoord, auxdata, stride, iwalker)
                ad_arr = numpy.asarray(ad_arr)

                # Test Functions... to be deleted...
                # ad_arr = numpy.array(list(iwalker.auxiliary_data[dataset])[:-2])[:,-1]
                # ad_arr = numpy.insert(ad_arr, 0, iwalker.pcoords[:,0][-1], axis=0)

                # pc_arr = iwalker.pcoords[:,1][-1]
                # closest_c = numpy.argmin(ad_arr)
                weight = iwalker.weight
                corr_assign = assign_file["statelabels"][
                    iwalker.iteration.summary.name - 1, iwalker.segment_summary.name
                ]
                if iwalker.iteration.summary.name == iteration_num:
                    # Dealing with cases where this is the first iteration we're looking at
                    # Taking only the first instance in tstate
                    term_frame_num = numpy.where(corr_assign == target_state_num)[0][0]
                    if trace_basis is False:
                        if source_state_num in corr_assign[: term_frame_num + 1]:
                            # Went from source to target in one iteration. neat.
                            source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                            indv_trace.append([iteration_num, segment_num, corr_assign[-1], *ad_arr, weight])
                            break
                    # Just a normal iteration where we reached target state and that's it.
                    indv_trace.append([iteration_num, segment_num, corr_assign[-1], *ad_arr, weight])

                elif iwalker.iteration.summary.name != iteration_num:
                    # Dealing with cases where we're in other iterations
                    if source_state_num not in corr_assign:
                        # The segment did not visit the source this iteration.
                        if iwalker.iteration.summary.name != iteration_num:
                            if target_state_num in corr_assign:
                                # Hey, this is a target > target transition.
                                # Breaking out...
                                run.close()
                                return None, None, None
                            else:
                                # If traj hasn't been in state A, and not in target iteration...
                                # add the whole iteration into list.
                                # Also making sure it's not a target -> target transition
                                indv_trace.append(
                                    [iwalker.iteration.summary.name, iwalker.segment_summary.name, corr_assign[-1],
                                     *ad_arr, weight]
                                )
                    else:
                        # This else captures cases where it was in the source in this iteration
                        # Looking for the last frame it was in source state
                        source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                        if target_state_num in corr_assign[source_frame_num:]:
                            # Catching the ouchie case where it went from source -> target -> target. Woopsies!
                            run.close()
                            return None, None, None
                        else:
                            # Final case where it's definitely source -> target
                            indv_trace.append(
                                [iwalker.iteration.summary.name, iwalker.segment_summary.name, corr_assign[-1],
                                 *ad_arr, weight]
                            )
                            if trace_basis is False:
                                break

            frame_info = []
            try:
                source_frame_num
            except NameError:
                source_frame_num = 0

            # Total number of frames necessary
            start_trace = (traj_len - source_frame_num) + ((len(indv_trace) - 1) * traj_len)
            end_trace = traj_len - term_frame_num
            frame_info.append(term_frame_num)
            for i in range(len(indv_trace) - 1, 1, -1):
                frame_info.append(traj_len)
            frame_info.append(source_frame_num)

            # Block for outputting the traj
            if out_traj:
                tqdm_bar.set_description(f"outputting traj for {iteration_num}.{segment_num}")
                if hdf5 is False:
                    trajectory = wa.BasicMDTrajectory(
                        traj_ext=out_traj_ext, state_ext=out_state_ext, top=out_top
                    )
                elif hdf5 is True:
                    trajectory = wa.HDF5MDTrajectory()
                else:
                    log.warning("unable to output trajectory")
                    return indv_trace, None, frame_info

                indv_traj = trajectory(trace)
                if trace_basis is False:
                    indv_traj = indv_traj[-start_trace:-end_trace]
            else:
                indv_traj = None

            # print(indv_trace)

            run.close()
            return indv_trace, indv_traj, frame_info

    def retain_succ(
            west_name=arguments.west_name,
            assign_name=arguments.assign_name,
            source_state_num=arguments.source_state_num,
            target_state_num=arguments.target_state_num,
            first_iter=arguments.first_iter,
            last_iter=arguments.last_iter,
            trace_basis=arguments.trace_basis,
            out_traj=arguments.out_traj,
            out_traj_ext=arguments.out_traj_ext,
            out_state_ext=arguments.out_state_ext,
            out_top=arguments.out_state_ext,
            out_dir=arguments.out_dir,
            hdf5=arguments.hdf5,
            rewrite_weights=arguments.rewrite_weights,
            use_ray=arguments.use_ray,
            threads=arguments.threads,
            pcoord=arguments.pcoord,
            auxdata=arguments.auxdata,
            output_name=arguments.extract_output,
            stride=arguments.stride,
    ):
        """
        Code that goes through an assign file (assign_name) and extracts iteration
        and segment number + its trace into a pickle object. Defaults to just
        the whole passage but can be overridden by trace_basis=False
        to just the transition time (between it last exits source to when
        it first touches target). Can also extract the trajectories along
        the way with out_traj=True.

        Arguments
        =========
        west_name : str, default: "west.h5"
            Name of the HDF5 File.

        assign_name : str, default: "assign.h5"
            Name of the `w_assign` output file.

        source_state_num : int, default: 0
            Index of the source state. Should match the state defined in
            ``west.cfg`` before running ``w_assign``.

        target_state_num : int, default: 1
            Index of the target state. Should match the state defined in
            ``west.cfg`` before running ``w_assign``.

        first_iter : int, default: 1
            First iteration to analyze. Inclusive.

        last_iter : int, default: 0
            Last iteration to analyze. Inclusive. 0 implies it will
            analyze all labeled iterations.

        trace_basis : bool, default: True
            Option to analyze each successful trajectory up till its
            basis state. False will only analyze the transition time
            (i.e. last it exited source until it first touches the target state).

        out_traj : bool, default: False
            Option to output trajectory files into ``out_dir``.

        out_traj_ext : str, default: ".nc"
            Extension of the segment files. The name of the file is assumed to be
            ``seg``, meaning the default name of the file is ``seg.nc``.

        out_state_ext : str, default: ".ncrst"
            Extension of the restart files. The name of the file is assumed to be
            ``seg``, meaning the default name the file is ``seg_nowat.ncrst``.

        out_top : str, default: "system.prmtop"
            Name of the topology file.
            Name is relative to ``$PWD/common_files/``.

        out_dir : str, default: "succ_traj"
            Name of directory to output the trajectories.
            Name is relative to ``$PWD``.

        hdf5 : bool, default: False
            Option to use the `HDF5MDTrajectory()` object in `westpa.analysis`
            instead. To be used with trajectories saved with the HDF5 Framework.

        rewrite_weights : bool, default: False
            Option to zero out the weights of all segments that are not part of
            the successful trajectory ensemble. Note this generates a new h5
            file with the ``_succ`` suffix added. Default name is thus ``west_succ.h5``.

        use_ray : bool, default: True if ray exists
            Option to turn ray on. This is assumed to be True if ray could be
            imported.

        threads : int, default: None
            Number of actors/workers for ray to initialize. It will take over all
            CPUs if None.

        pcoord : bool, default: True
            Boolean confirming if you would like to include the progress coordinate
            of the last frame in the ``output.pickle`` file. Default to True.

        auxdata : list of strings, default: None
            Auxiliary data set you would like to include in the ``output.pickle`` file.
            None means you don't want any. Only includes the last frame.

        output_name : str, default: output.pickle
            Name of the output pickle file.

        stride : int, default: 1
            The number of frames to output per segment.

        Notes
        =====
        The following files are saved/outputted to disk.

        'output.pickle' : pickle obj
            A list of the form [n_traj, n_frame, [n_iter, n_seg, state_id, [pcoord/auxdata], weight]]. Note that each
            list runs backwards from the target iteration.

        'frame_info.pickle': pickle obj
            A list of the form [n_traj, [index number of each segment]]. Note that
            each list runs backwards from the target iteration.
        """
        # Variables validation
        if last_iter == 0:
            with h5py.File(assign_name, "r") as assign_file:
                last_iter = len(assign_file["nsegs"])
        else:
            assert type(last_iter) is int, "last_iter is not legal and must be int."

        # Create Output file
        try:
            mkdir(out_dir)
        except FileExistsError:
            print(f"Folder {out_dir} already exists. Files within might be overwritten.")

        # Copying the file
        name_root = west_name.rsplit(".h5", maxsplit=1)[0]
        new_file = f"{out_dir}/{name_root}_succ.h5"
        # if not exists(new_file): # Always recopy file...
        copyfile(west_name, new_file)

        # Prepping final list to be outputted
        trace_out_list = []
        frame_info_list = []

        # Ray stuff...
        if use_ray:
            ray.init()
            if threads == 0:
                n_actors = int(ray.available_resources().get("CPU", 1))
            else:
                n_actors = threads
            all_ray_actors = []
            for i in range(n_actors):
                all_ray_actors.append(TraceActor.remote(assign_name, new_file))

            # Yes, tracing backwards from the last iteration. This will (theoretically) allow us to catch
            # duplicates more efficiently.
            with h5py.File(assign_name, "r") as assign_file:
                tqdm_iter = trange(last_iter, first_iter - 1, -1, desc="iter")
                for n_iter in tqdm_iter:
                    all_ray_tasks = []
                    for n_seg in range(assign_file["nsegs"][n_iter - 1]):
                        if target_state_num in set(assign_file["statelabels"][n_iter - 1, n_seg]):
                            all_ray_tasks.append(all_ray_actors[n_seg % n_actors].trace_seg_to_last_state.remote(
                                source_state_num,
                                target_state_num,
                                n_iter,
                                n_seg,
                                first_iter,
                                trace_basis,
                                out_traj,
                                out_traj_ext,
                                out_state_ext,
                                out_top,
                                hdf5,
                                pcoord,
                                auxdata,
                                stride,
                            ))

                    while all_ray_tasks:
                        finished, all_ray_tasks = ray.wait(all_ray_tasks, num_returns=min(n_actors * 5,
                                                                                          len(all_ray_tasks)))
                        results = ray.get(finished)
                        for each_return in results:
                            (trace_output, traj_output, frame_info) = each_return
                            if trace_output is not None:
                                trace_out_list.append(trace_output)
                            if traj_output is not None:
                                traj_output.save(f"{out_dir}/{n_iter}_{n_seg}{out_traj_ext}")
                            if frame_info is not None:
                                frame_info_list.append(frame_info)

            # Shutting down ray since we're done with parallelization
            ray.shutdown()
        else:
            # Yes, tracing backwards from the last iteration. This will (theoretically) allow us to catch
            # duplicates more efficiently.
            with h5py.File(assign_name, "r") as assign_file:
                tqdm_iter = trange(last_iter, first_iter - 1, -1, desc="iter")
                for n_iter in tqdm_iter:
                    for n_seg in range(assign_file["nsegs"][n_iter - 1]):
                        if target_state_num in set(assign_file["statelabels"][n_iter - 1, n_seg]):
                            trace_output, traj_output, frame_info = trace_seg_to_last_state(
                                source_state_num,
                                target_state_num,
                                new_file,
                                assign_file,
                                n_iter,
                                n_seg,
                                trace_basis,
                                out_traj,
                                out_traj_ext,
                                out_state_ext,
                                out_top,
                                hdf5,
                                tqdm_iter,
                                pcoord,
                                auxdata,
                                stride,
                            )

                            if trace_output is not None:
                                trace_out_list.append(trace_output)
                            if traj_output is not None:
                                traj_output.save(f"{out_dir}/{n_iter}_{n_seg}{out_traj_ext}")
                            if frame_info is not None:
                                frame_info_list.append(frame_info)

        # Output list
        trace_out_list = sorted(trace_out_list, key=lambda x: (-x[0][0], x[0][1]))

        with open(f"{out_dir}/{output_name}", "wb") as fo:
            pickle.dump(trace_out_list, fo)

        with open(f"{out_dir}/frame_info.pickle", "wb") as fo:
            pickle.dump(frame_info_list, fo)

        # Finally, zero out (iter,seg) that do not fall in this "successful" list.
        if rewrite_weights:
            exclusive_set = {tuple([pair[0], pair[1]]) for ilist in trace_out_list for pair in ilist}
            with h5py.File(new_file, "r+") as h5file, h5py.File(assign_name, "r") as assign_file:
                for n_iter in tqdm_iter:
                    for n_seg in range(assign_file["nsegs"][n_iter - 1]):
                        if (n_iter, n_seg) not in exclusive_set:
                            h5file[f"iterations/iter_{n_iter:>08}/seg_index"]["weight", n_seg] = 0

    retain_succ()


def main(arguments):
    """
    Main function that executes the `match` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    if arguments.we:
        log.debug('Running Extract WE mode.')
        we(arguments)
    else:
        log.debug('Running Extract Standard mode.')
        standard(arguments)


if __name__ == "__main__":
    """
    For calling `extract.py` directly. Note all of the parameters are specified manually here.
    """
    # TODO: This is broken at the moment.
    import argparse

    args = argparse.Namespace(
        we=True,  # Analyzing a WE simulation
        west_name="west.h5",  # Name of input HDF5 file (e.g., west.h5)
        assign_name="ANALYSIS/C7_EQ/assign.h5",  # Name of input assign.h5 file
        source_state_num=0,  # Index of the source state as defined in discretize.
        target_state_num=1,  # Index of the target state as defined in discretize.
        first_iter=1,  # First iteration to analyze. Inclusive
        last_iter=200,  # Last iteration to analyze. Inclusive. 0 implies it will analyze all labeled iterations.
        trace_basis=True,  # Option to analyze each successful trajectory up till its basis state.
        out_traj=False,  # Option to output trajectory files into `out_dir`. Will take much longer.
        out_traj_ext=".nc",  # Extension of the segment files. Defaults to `seg{out_traj_ext}`.
        out_state_ext=".ncrst",  # Extension of the restart files. Defaults to `seg{out_state_ext}`.
        out_top="system.prmtop",  # Name of the parameter file. Name relative to `$WEST_SIM_ROOT/common_files`.
        out_dir="succ_traj",  # Name of directory to output the trajectories.
        extract_output="output.pickle",  # Name of the pickle file to be outputted.
        hdf5=False,  # Enable if trajectories are saved with the HDF5 Framework in WESTPA.
        rewrite_weights=False,  # Option to zero out the weights of all segments that are not a successful trajectory.
        pcoord=True,  # Option to output the pcoord into the `output.pickle`.
        auxdata=["phi", "psi"],  # Additional auxiliary data to save into `output.pickle`.
        use_ray=True,  # Enable Ray.
        threads=0,  # How many Ray threads/actors to use. Defaults to 0, which wil use all auto-detected resources.
        stride=1,  # How much to stride the dataset
    )
    main(args)
