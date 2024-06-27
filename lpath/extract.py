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
# Defaults to the barrier crossing time (from last exit from source to first
#     touch target) for all of the above.
# Make ``trace_basis`` True to trace only the barrier crossing time.

import pickle
from collections import Counter
from copy import deepcopy
from shutil import copyfile

import numpy
from tqdm.auto import tqdm, trange

from lpath.io import load_file, expanded_load, EmptyOutputError
from lpath.extloader import *

from ._logger import Logger

log = Logger().get_logger(__name__)


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
        The closest index (in indices) to ref_value.

    """
    return min([v for v in indices if v >= ref_value] or [None])


def assign_color_frame(source_indices, target_indices):
    """
    Assign color to each target frame.

    Parameters
    ----------
    source_indices : list or numpy.ndarray
        A list or array of indices assigned source state.

    target_indices : list or numpy.ndarray
        A list or array of indices assigned sink state.

    Returns
    -------
    target_colors : dictionary
        A list mapping each target index in ``target_indices`` with a frame it should trace back to.
        Essentially the frame when the color of the trajectory switched from sink to source.

    """
    # Assign default color to each target index. 0 by default.
    target_colors = dict()
    last_source = 0
    test_source = 0
    for idx, val in enumerate(target_indices):
        if val > test_source:
            last_source = test_source
        target_colors[val] = last_source
        # Finding the next source after each target, which should be the frame saved for the next idx.
        # At the last loop, last_source should be None.
        test_source = find_min_distance(val, source_indices)

    assert len(target_colors) == len(target_indices)

    return target_colors


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

    reduced_keys = [key for (key, count) in full_count.items() if count > 1]

    # For debugging purposes.
    # log.debug(f"Full Count: {full_count}")
    # log.debug(f"Reduced Keys: {reduced_keys}")
    # Running backwards so indices are maintained.
    for delete in tqdm(reduced_keys[::-1], desc='Cleaning redundant transitions', leave=False):
        # Determine indices of where the duplicates happen
        pop_list = numpy.argwhere(output_array[:, 1] == delete).flatten()
        pop_list.sort()
        log.debug(f"All indices where {delete} occur: {pop_list}")
        # Remove them except for the last instance
        for j in pop_list[::-1][1:]:
            input_array.pop(j)

    output_array = numpy.asarray(input_array)
    log.debug(f'Indices of all true successful transitions: {output_array}')

    return output_array


def count_tmatrix_row(source_index, trajectory, n_states, source_num, target_num):
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
            if jstate in (source_num, target_num):
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
    #     f"None of your target state assignments {target_index} occur after the "
    #      "last visit to the source state {source_index}."

    # Now do the calculations
    transitions = []
    for idx in tqdm(source_indices, desc='Tracing successful transitions', leave=False):
        check = find_min_distance(idx, target_indices)
        if check:
            transitions.append([idx, check])

    log.debug(f'Indices of all potentially successful transitions: {transitions}')

    return source_indices, target_indices, transitions


def raise_warnings(output_array, statistics):
    """
    Raise warnings and Errors towards common failure modes.

    Parameters
    ----------
    output_array : list
        A list of lists containing traced trajectories

    statistics : bool, default: False
        A flag to report statistics.

    """
    # If no successful trajectories...
    if len(output_array) == 0:
        raise EmptyOutputError

    lengths = [len(a) for a in output_array]
    if min(lengths) < 10:
        log.warning('Extracted trajectories have a minimum length of <10 frames, which will '
                    'affect the quality of pattern matching further downstream. For standard simulations, use a '
                    'smaller ``stride`` to retain more data points. For WE simulations, use a larger ``stride`` to '
                    'output more sub-tau values. Another option is to use ``--trace-basis`` to include trajectory '
                    'history all the way back to the "basis state". Proceeding to output pickle object anyway.')

    sets = set([i[2] for a in output_array for i in a])
    log.info(f"States {sets} visited within all successful trajectories.")
    if len(sets) < 5:
        log.warning(f'Your trajectories visited only {len(sets)}')

    if statistics:
        log.info(f'Mean length ± 1σ of trajectory: {numpy.average(lengths)} ± {numpy.std(lengths)}\n')


def create_pickle_obj(transitions, states, weight, features=None):
    """
    Main function that transforms a list of frame transitions into the pickle object. For standard simulations only.

    Parameters
    ----------
    transitions : list or numpy.ndarray
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
        ad_arr = [[] for _ in range(len(states))]
    else:
        if isinstance(features, numpy.ndarray):
            ad_arr = features.tolist()
        elif isinstance(features, list):
            ad_arr = features
        else:
            ad_arr = list(features)
        # Assuming there are at least 2 frames, here.
        # We want shape of [[1],[2],[3],...] not [[1,2,3,...]].
        try:
            isinstance(ad_arr[1], list)
        except IndexError:
            ad_arr = [[i] for i in ad_arr[0]]

    output_list = []

    # Going through each transition
    for idx, transition in enumerate(transitions):
        indv_trace = []
        # Add all the frames between target + source. This is controlled by stride, which dictates how often you load.
        for j in range(int(transition[0]), int(transition[1] + 1)):
            indv_trace.append([1, idx, states[j], *ad_arr[j], j, weight])
        output_list.append(deepcopy(indv_trace))

    return output_list


def standard(arguments):
    """
    Main function that executes the whole standard ``extract`` procedure.
    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    input_array = load_file(arguments.extract_input, arguments.stride)
    n_states = len(input_array) - 1

    if arguments.pcoord is True:
        if arguments.featurization_name is None:
            log.warning('The --pcoord flag is specified but no file is '
                        'specified with --extract-pcoord. Skipping output.')
            features = None
        else:
            features = [expanded_load(arguments.featurization_name, arguments.feature_stride)]
    else:
        features = None

    source_index, target_index, transitions = find_transitions(input_array, arguments.source_state_num,
                                                               arguments.target_state_num)

    new_transitions = clean_self_to_self(transitions)

    if arguments.trace_basis:
        target_color_dict = assign_color_frame(source_index, target_index)
        for transition in new_transitions:
            transition[0] = target_color_dict[transition[1]]

    weight = count_tmatrix_row(source_index, input_array, n_states, arguments.source_state_num,
                               arguments.target_state_num)

    # Generate and write pickle object.
    final_obj = create_pickle_obj(new_transitions, input_array, weight / len(transitions), features)

    raise_warnings(final_obj, arguments.stats)

    with open(arguments.extract_output, "wb") as fo:
        pickle.dump(final_obj, fo)


# Here are functions for WE.
def we(arguments):
    """
    Main function that executes the whole WE ``extract`` procedure.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    import h5py
    import westpa.analysis as wa

    def process_ad_arr(pcoord, auxdata, iwalker):
        """
        Inner function that deals with preparing ad_arr

        Parameters
        ----------
        pcoord : bool
            Marks whether to add progress coordinate in ad_arr or not.

        auxdata : list or None
            A list of auxiliary datasets to add or None

        iwalker : westpa.analysis.core.Walker
            A walker object from ``westpa.analysis`` of relevant segment.

        Returns
        -------
        ad_arr : list
            A list of relevant data to be outputted.

        """
        total_frames = iwalker.pcoords.shape[0]
        if pcoord is True:
            ad_arr = iwalker.pcoords.tolist()
        else:
            ad_arr = []

        if auxdata is not None:
            if len(auxdata) == 0:
                # auxdata is set to 0 when auxall is called, so grabbing all the aux datasets.
                auxdata = list(iwalker.auxiliary_data.keys())
                log.info(f'Exporting the following datasets: {auxdata}')
            for dataset_name in auxdata:
                # Using list comprehension, since numpy.append is weird at times.
                if len(iwalker.auxiliary_data[dataset_name].shape) > 1:
                    _ = [row.append(item) for (row, dataset) in zip(ad_arr, iwalker.auxiliary_data[dataset_name])
                         for item in dataset]
                else:
                    _ = [row.append(dataset) for (row, dataset) in zip(ad_arr, iwalker.auxiliary_data[dataset_name])]

        if len(ad_arr) == 0:
            ad_arr = [[] for _ in range(total_frames)]

        # Turning it back to numpy array for easier indexing.
        ad_arr = numpy.asarray(ad_arr, dtype=object)

        return ad_arr

    def frame_range(source_frame_num, term_frame_num, n_frames, stride_step):
        """
        Generator the frames we need to output based on ``stride_step`` and ``trace_basis``.
        Basically a replacement ``range`` object with extra checks for outputting.


        Parameters
        ----------
        source_frame_num : int
            Index of where trajectory last exited source state.

        term_frame_num : int
            Index of where trajectory first reached target state.

        n_frames : int
            Number of frames in each segment

        stride_step : int
            How often to output frames

        Yields
        ------
        frame_loop : list or range
            frames to output, outputted in reverse order as list or range object.

        """
        base_frames = range(n_frames - 1, -1, -stride_step)

        if (source_frame_num, term_frame_num) != (-1, n_frames):
            base_frames = list(base_frames)
            if source_frame_num != -1:
                # Source frame and every base frame after
                if source_frame_num not in base_frames:
                    base_frames = base_frames + [source_frame_num]
                base_frames = sorted([frame for frame in base_frames if frame >= source_frame_num], reverse=True)

            if term_frame_num != n_frames:
                # Every base frame before (and including) term frame
                if term_frame_num not in base_frames:
                    base_frames = base_frames + [term_frame_num]
                base_frames = sorted([frame for frame in base_frames if frame <= term_frame_num], reverse=True)

        return base_frames

    def check_west_assign(west_name, assign_name, first_iter, last_iter):
        """
        A function to check if the assign.h5 corresponds with west.h5.

        Parameters
        ----------
        west_name : str
            Name of ``west.h5`` file. Identical from what's inputted to the argparser.

        assign_name : str
            Name of ``assign.h5`` file. Identical from what's inputted to the argparser.

        first_iter : int or float
            First iteration number. 1-based.

        last_iter : int or float
            Last iteration number. 1-based.

        """
        with h5py.File(west_name) as west_file, h5py.File(assign_name) as assign_file:
            for iter_num in range(first_iter, last_iter):
                assert west_file['summary']['n_particles', iter_num] == assign_file['nsegs'][iter_num], (
                    f"The assign.h5 file provided doesn't match the west.h5 file provided at iteration {iter_num +1}.")

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
                    # Grabbing relevant datasets
                    ad_arr = process_ad_arr(pcoord, auxdata, iwalker)
                    weight = iwalker.weight
                    corr_assign = assign_file["statelabels"][
                        iwalker.iteration.summary.name - 1, iwalker.segment_summary.name
                    ]

                    # Stride parameters
                    total_frames = iwalker.pcoords.shape[0]
                    stride_step = total_frames // stride or 1

                    if iwalker.iteration.summary.name == iteration_num:
                        # Dealing with cases where this is the first iteration we're looking at
                        # If multiple, taking only the first instance we reached target state
                        term_frame_num = numpy.where(corr_assign == target_state_num)[0][0]
                        if source_state_num in corr_assign[:term_frame_num + 1]:
                            # Went from source to target in one iteration. neat.
                            # Grabbing last instance it exited source state
                            source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                            # print(source_frame_num, term_frame_num)
                            frame_loop = frame_range(source_frame_num, term_frame_num, total_frames, stride_step)
                            # print(frame_loop)
                            for frame_index in frame_loop:
                                indv_trace.append([iteration_num, segment_num, corr_assign[frame_index],
                                                   *ad_arr[frame_index], frame_index, weight])
                            if trace_basis is False:
                                break
                        else:
                            # Just a normal iteration where we reached target state. Output everything in stride.
                            frame_loop = frame_range(-1, term_frame_num, total_frames, stride_step)
                            for frame_index in frame_loop:
                                indv_trace.append([iteration_num, segment_num, corr_assign[frame_index],
                                                   *ad_arr[frame_index], frame_index, weight])
                    elif iwalker.iteration.summary.name != iteration_num:
                        # Dealing with cases where we're in other iterations
                        if source_state_num not in corr_assign:
                            # The segment did not visit the source this iteration.
                            if iwalker.iteration.summary.name != iteration_num:
                                if target_state_num in corr_assign:
                                    # Hey, this is a target > target transition.
                                    # Breaking out...
                                    return None, None
                                else:
                                    # If traj hasn't been in source state, and not in target iteration...
                                    # add the whole iteration into list.
                                    # Also making sure it's not a target -> target transition
                                    for frame_index in frame_range(-1, total_frames, total_frames, stride_step):
                                        indv_trace.append(
                                            [iwalker.iteration.summary.name, iwalker.segment_summary.name,
                                             corr_assign[frame_index], *ad_arr[frame_index], frame_index, weight]
                                        )
                        else:
                            # This else captures cases where it was in the source in this iteration
                            # Looking for the last frame it was in source state
                            source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                            if target_state_num in corr_assign[source_frame_num:]:
                                # Catching the ouchie case where it went from source -> target -> target. Whoopsies!
                                return None, None
                            else:
                                # Final case where it's definitely source -> target
                                frame_loop = frame_range(source_frame_num, total_frames, total_frames, stride_step)
                                for frame_index in frame_loop:
                                    indv_trace.append(
                                        [iwalker.iteration.summary.name, iwalker.segment_summary.name,
                                         corr_assign[frame_index], *ad_arr[frame_index], frame_index, weight]
                                    )
                                if trace_basis is False:
                                    break

                try:
                    source_frame_num
                except NameError:
                    source_frame_num = 0
                # Total number of frames necessary
                start_trace = (traj_len - source_frame_num) + ((len(indv_trace) - 1) * traj_len)
                end_trace = traj_len - term_frame_num

                # Block for outputting the traj
                if out_traj:
                    if hdf5 is False:
                        trajectory = wa.BasicMDTrajectory(
                            traj_ext=out_traj_ext, state_ext=out_state_ext, top=out_top
                        )
                    elif hdf5 is True:
                        trajectory = wa.HDF5MDTrajectory()
                    else:
                        log.info("Unable to output trajectory")
                        return indv_trace, None

                    # Extra check such that it won't output past first_iter.
                    if first_iter != 1:
                        max_length = iteration_num - first_iter + 1
                        trace = run.iteration(iteration_num).walker(segment_num).trace(max_length=max_length)

                    indv_traj = trajectory(trace)
                    if trace_basis is False:
                        indv_traj = indv_traj[-start_trace:-end_trace]
                else:
                    indv_traj = None

                # print(f'{indv_trace}')

                return indv_trace, indv_traj
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
                # Grabbing relevant datasets
                ad_arr = process_ad_arr(pcoord, auxdata, iwalker)
                weight = iwalker.weight
                corr_assign = assign_file["statelabels"][
                    iwalker.iteration.summary.name - 1, iwalker.segment_summary.name
                ]

                # Stride parameters
                total_frames = iwalker.pcoords.shape[0]
                stride_step = total_frames // stride or 1

                if iwalker.iteration.summary.name == iteration_num:
                    # Dealing with cases where this is the first iteration we're looking at
                    # If multiple, taking only the first instance we reached target state
                    term_frame_num = numpy.where(corr_assign == target_state_num)[0][0]
                    if source_state_num in corr_assign[:term_frame_num + 1]:
                        # Went from source to target in one iteration. neat.
                        # Grabbing last instance it exited source state
                        source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                        frame_loop = frame_range(source_frame_num, term_frame_num, total_frames, stride_step)
                        for frame_index in frame_loop:
                            indv_trace.append([iteration_num, segment_num, corr_assign[frame_index],
                                               *ad_arr[frame_index], frame_index, weight])
                        if trace_basis is False:
                            break
                    else:
                        # Just a normal iteration where we reached target state. Output everything in stride.
                        frame_loop = frame_range(-1, term_frame_num, total_frames, stride_step)
                        for frame_index in frame_loop:
                            indv_trace.append([iteration_num, segment_num, corr_assign[frame_index],
                                               *ad_arr[frame_index], frame_index, weight])
                elif iwalker.iteration.summary.name != iteration_num:
                    # Dealing with cases where we're in other iterations
                    if source_state_num not in corr_assign:
                        # The segment did not visit the source this iteration.
                        if iwalker.iteration.summary.name != iteration_num:
                            if target_state_num in corr_assign:
                                # Hey, this is a target > target transition.
                                # Breaking out...
                                return None, None
                            else:
                                # If traj hasn't been in state A, and not in target iteration...
                                # add the whole iteration into list.
                                # Also making sure it's not a target -> target transition
                                for frame_index in frame_range(-1, total_frames, total_frames, stride_step):
                                    indv_trace.append(
                                        [iwalker.iteration.summary.name, iwalker.segment_summary.name,
                                         corr_assign[frame_index], *ad_arr[frame_index], frame_index, weight]
                                    )
                    else:
                        # This else captures cases where it was in the source in this iteration
                        # Looking for the last frame it was in source state
                        source_frame_num = numpy.where(corr_assign == source_state_num)[0][-1]
                        if target_state_num in corr_assign[source_frame_num:]:
                            # Catching the ouchie case where it went from source -> target -> target. Whoopsies!
                            return None, None
                        else:
                            # Final case where it's definitely source -> target
                            frame_loop = frame_range(source_frame_num, total_frames, total_frames, stride_step)
                            for frame_index in frame_loop:
                                indv_trace.append(
                                    [iwalker.iteration.summary.name, iwalker.segment_summary.name,
                                     corr_assign[frame_index], *ad_arr[frame_index], frame_index, weight]
                                )
                            if trace_basis is False:
                                break

            try:
                source_frame_num
            except NameError:
                source_frame_num = 0
            # Total number of frames necessary
            start_trace = (traj_len - source_frame_num) + ((len(indv_trace) - 1) * traj_len)
            end_trace = traj_len - term_frame_num

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
                    return indv_trace, None

                indv_traj = trajectory(trace)
                if trace_basis is False:
                    indv_traj = indv_traj[-start_trace:-end_trace]
            else:
                indv_traj = None

            # print(indv_trace)

            run.close()
            return indv_trace, indv_traj

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
            exclude_short=arguments.exclude_short,
    ):
        """
        Code that goes through an assign file (assign_name) and extracts iteration
        and segment number + its trace into a pickle object. Defaults to just
        the transition time (between it last exits source to when
        it first touches target) but can be switched to the whole passage (all the way
        to the basis state) with ``trace_basis=True``. Can also extract the full
        trajectories along the way with out_traj=True.

        Arguments
        =========
        west_name : str, default: "west.h5"
            Name of the HDF5 File.

        assign_name : str, default: "assign.h5"
            Name of the ``w_assign`` output file.

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

        trace_basis : bool, default: False
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
            Option to use the ``HDF5MDTrajectory()`` object in ``westpa.analysis``
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

        output_name : str, default: 'succ_traj/output.pickle'
            Name of the output pickle file.

        stride : int, default: 1
            The number of frames to output per segment.

        Notes
        =====
        The following files are saved/outputted to disk.

        'output.pickle' : pickle obj
            A list of the form [n_traj, n_frame, [n_iter, n_seg, state_id, [pcoord/auxdata], weight]]. Note that each
            list runs backwards from the target iteration.

        """
        # Copying the file
        name_root = west_name.rsplit(".h5", maxsplit=1)[0]
        new_file = f"{out_dir}/{name_root}_succ.h5"
        # if not exists(new_file): # Always recopy file...
        copyfile(west_name, new_file)

        # Prepping final list to be outputted
        trace_out_list = []

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
                tqdm_iter = trange(last_iter, first_iter - 1, -1, desc="iterations")
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
                            (trace_output, traj_output) = each_return
                            if trace_output:
                                trace_out_list.append(trace_output[::-1])
                            if traj_output:
                                traj_output.save(f"{out_dir}/{n_iter}_{n_seg}{out_traj_ext}")

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
                            trace_output, traj_output = trace_seg_to_last_state(
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

                            if trace_output:
                                trace_out_list.append(trace_output[::-1])
                            if traj_output:
                                traj_output.save(f"{out_dir}/{n_iter}_{n_seg}{out_traj_ext}")

        # Failure modes and warnings... and result stats!
        raise_warnings(trace_out_list, arguments.stats)

        # Output list
        trace_out_list = sorted(trace_out_list, key=lambda x: (-x[-1][0], x[-1][1]))

        if exclude_short:
            del_list = []
            for idx, pathway in enumerate(trace_out_list):
                if len(pathway) < exclude_short:
                    del_list.append(idx)

            if len(del_list) > 0:
                for jdx in del_list[::-1]:
                    # Deleting from larger to lower
                    deleted = trace_out_list.pop(jdx)
                    log.debug(f'Deleting pathway index {jdx}: {deleted}')
                log.info(f'Indices of pathway removed: {del_list}.')
                log.info(f'Removed {len(del_list)} trajectories of length < {exclude_short} frames.')

        with open(f"{output_name}", "wb") as fo:
            pickle.dump(trace_out_list, fo)

        # Finally, zero out (iter,seg) that do not fall in this "successful" list.
        if rewrite_weights:
            exclusive_set = {tuple([pair[0], pair[1]]) for ilist in trace_out_list for pair in ilist}
            with h5py.File(new_file, "r+") as h5file, h5py.File(assign_name, "r") as afile:
                for n_iter in tqdm_iter:
                    for n_seg in range(afile["nsegs"][n_iter - 1]):
                        if (n_iter, n_seg) not in exclusive_set:
                            h5file[f"iterations/iter_{n_iter:>08}/seg_index"]["weight", n_seg] = 0

    # Variables validation
    if arguments.last_iter == 0:
        with h5py.File(arguments.assign_name, "r") as assign_file:
            setattr(arguments, 'last_iter', len(assign_file["nsegs"])-1)
    else:
        assert type(arguments.last_iter) is int, "last_iter is not legal and must be int."

    check_west_assign(arguments.west_name, arguments.assign_name, arguments.first_iter, arguments.last_iter)
    retain_succ(last_iter=arguments.last_iter)


def main(arguments):
    """
    Main function that executes the ``match`` step.

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
