Usage
=====

This page details how to use ``lpath``.  Follow the `Getting Started`_ page to learn how to install ``lpath``.

.. _Getting Started: https://lpath.readthedocs.io/en/latest/getting_started.html


Introduction
------------
``lpath`` is made of four separate steps: ``discretize``, ``extract``, ``match``, and ``plot``. All these steps take a variety of different parameters, which are all listed in the `API`_ page.

1. The ``discretize`` step is used to discretize your trajectories for the purpose of finding successful transitions.
2. The ``extract`` step takes what's assigned in ``discretize`` and identifies all instances where there is a successful transition.
3. The ``match`` step takes what's outputted in ``extract`` and cross pattern match to identify pathway classes. It is possible to reassign states in this step.
4. The ``plot`` step takes what's outputted in ``match`` and allows you to easily make different plots. The LPATHPlot class objects contains many pre-calculated datasets for custom plotting.


There are two different ways of running these steps. Due to the sheer amount of parameter options, it is recommended that users start with the Jupyter notebook.

1. Import each step's ``main()`` function and run everything in an interactive python session (e.g., Jupyter notebook).  **[RECOMMENDED]**
2. Run through the command line (e.g., ``lpath discretize -I west.h5 --assign-arguments='--config-from-file --scheme TEST'``.


.. _API: https://lpath.readthedocs.io/en/latest/api.html

Molecular Dynamics
------------------
``lpath`` is written to cluster pathways from molecular dynamics simulation. The input file will be a numpy or text file of the features used to assign states. For example, an alanine dipeptide system will use the phi, psi angles.

Discretize
__________
In this step, we will assign each frame of an MD trajectory based on the phi/psi dihedral angles saved in ``dihedral.npy``.
Here is an example ``assign`` function (in a file called ``module.py``) used for assigning states::

    def assign_dih(input_array):
        """
        This is an example function for mapping a list of features to state IDs. This should be subclassed.

        Parameters
        ----------
        input_array : numpy.ndarray
            An array generated from load_file.

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


We will monkey-patch this function into ``lpath``.

1. From the command line, run the following::

    lpath discretize -I dihedral.npy -O states.npy -af module.assign_dih --stride 100


2. We've read in ``dihedral.npy`` with a stride step of 100. This smaller dataset will be discretized with module.assign_dih. This will generate a ``states.npy`` file to be used in the ``extract`` step.

Extract
_______
In this step, we will identify any successful transitions in the trajectory. We will be looking at the C7\ :sub:`eq` to C7\ :sub:`ax` transition.
Since we already read in the data every 100 frames with ``--stride 100`` in the ``discretize`` step, we `do not` need to use ``--stride`` again.

Users who haven't used stride in ``discretize`` and would like their dataset reduced at this point may use ``--stride`` here.

Note: If you read in your trajectory in ``discretize`` with a different stride, your other pre-calculated datasets might no longer match frame-by-frame. For that dataset to be used for reassignment later on (i.e. in ``match``), use the ``--feature-stride 100`` option and provide the appropriate stride step.

1. From the command line, run the following::

    lpath extract --extract-input states.npy --extract-output pathways.pickle --source-state 1 --target-state 2


Match
_____
In this step, we will pattern match any successful transitions we've identified in ``extract``. We will, again, be looking at the C7\ :sub:`eq` to C7\ :sub:`ax` transition.

1. From the command line, run the following::

    lpath match --input-pickle succ_traj/pathways.pickle --output-pickle succ_traj/match-output.pickle \
    --cluster-labels-output succ_traj/cluster_labels.npy

2. After the comparison process is completed, it should show you the dendrogram. Closing the figure should trigger prompts to guide you further.

3. Input ``y`` if you think the threshold (horizontal line which dictates how many clusters there are) should be at a different value. Otherwise, input ``n`` and tell the program how many clusters you want at the end.

Plot
____
This step will help you plot some of the most common graphs, such as dendrograms and histograms, directly from the pickle object generated from match. Users may also elect to use the plotting scripts from the ``examples`` folder.
There is a script to plot ``NetworkX`` plots there.

More specifically, the following graphs will be made in the ``plots`` folder::

* Dendrogram showing separation between clusters
* Weights/Cluster bar graph
* Target iteration histograms (per cluster)
* Event duration histograms (per cluster)


From the command line, run the following and it should generate a separate file for each of the above graphs::

    lpath plot --plot-input succ_traj/match-output.pickle

More options for customizing the graphs can be found by running ``lpath plot --help``.

Weighted Ensemble Simulations
-----------------------------
``lpath`` is written to cluster pathways generated by the `WESTPA`_ software suite. Make sure `WESTPA`_ is installed. See the `Getting Started`_ page for more information.

.. _WESTPA: https://westpa.github.io/

Discretize
__________
We will use `WESTPA`_'s ``w_assign`` tool to assign to states. See the tool's `wiki`_ page and `Sphinx`_ documentation for more information about the tool.

.. _wiki: https://github.com/westpa/westpa/wiki/man:w_assign
.. _Sphinx: https://westpa.readthedocs.io/en/latest/documentation/cli/w_assign.html


We'll try to discretize a ``multi.h5`` (generated with ``w_multi_west --ibstates``) with ``w_assign`` based on what's defined with the ``TEST`` scheme in the ``west.cfg``. In ``TEST``, a rectilinear grid is constructed in the progress coordinate space and certain "bins" are selected and assigned to each state.

1. Run the following in the command line to run ``w_assign``::

    lpath discretize -we -W multi.h5 -A ANALYSIS/TEST/assign.h5 \
        --assign-args="-W multi.h5 -r west.cfg --config-from-file --scheme TEST"


Extract
_______
In this step, we will identify any successful transitions in the trajectory. We will be looking at the C7\ :sub:`eq` to C7\ :sub:`ax` transition.
If you are looking to compare using segment IDs in the next step (not recommended for simulations combined with ``w_multi_west``) or want to include the waiting time (time spent in the source state) in the pattern matching, make sure you turn on ``--trace-basis`` to trace all the way back to the basis state. Do note that this significantly increases the time it requires to extract all successful trajectories.

1. From the command line, run the following::

    lpath extract -we -W multi.h5 -A ANALYSIS/TEST/assign.h5 --source-state 1 \
        --target-state 2 --extract-output output.pickle --out-dir succ_traj


Match
_____
In this step, we will pattern match any successful transitions we've identified in ``extract``. We will, again, be looking at the C7\ :sub:`eq` to C7\ :sub:`ax` transition.
This will do the pattern matching and output individual h5 files for each cluster.

1. From the command line, run the following::

    lpath match -we --input-pickle succ_traj/output.pickle --output-pickle succ_traj/match-output.pickle  --cluster-labels-output succ_traj/cluster_labels.npy \
        --export-h5 --file-pattern "west_succ_c{}.h5"

2. After the comparison process is completed, it should show you the dendrogram. Closing the figure should trigger prompts to guide you further.

3. Input ``y`` if you think the threshold (horizontal line which dictates how many clusters there are) should be at a different value. Otherwise, input ``n`` and tell the program how many clusters you want at the end.


For cases where you want to run pattern matching comparison between segment IDs, you will have to use the largest common substring ``--substring`` option. By default, the longest common subsequence algorithm is used.::

    lpath match -we --input-pickle succ_traj/output.pickle --output-pickle succ_traj/match-output.pickle --cluster-labels-output succ_traj/cluster_labels.npy \
        --export-h5 --file-pattern "west_succ_c{}.h5" --reassign-method "reassign_segid" --substring


Plot
____
This step will help you plot some of the most common graphs, such as dendrograms and histograms, directly from the pickle object generated from match. Users may also elect to use the plotting scripts from the ``examples`` folder.
There is a script to plot ``NetworkX`` plots there.

More specifically, the following graphs will be made in the ``plots`` folder::

* Dendrogram showing separation between clusters
* Weights/Cluster bar graph
* Target iteration histograms (per cluster)
* Event duration histograms (per cluster)


From the command line, run the following and it should generate a separate file for each of the above graphs::

    lpath plot --plot-input succ_traj/match-output.pickle

More options for customizing the graphs can be found by running ``lpath plot --help``.


Example Reassign file
---------------------

The following is a reassign function if you decides to reclassify your states::

    def reassign_custom(data, pathways, dictionary, assign_file=None):
        """
        Reclassify/assign frames into different states. This is highly
        specific to the system. If w_assign's definition is sufficient,
        you can proceed with what's made in the previous step
        using ``reassign_identity``.

        In this example, the dictionary maps state idx to its corresponding ``state_string``.
        We suggest using alphabets as states.

        Parameters
        ----------
        data : list
            An array with the data necessary to reassign, as extracted from ``output.pickle``.

        pathways : numpy.ndarray
            An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

        dictionary : dict
            An empty dictionary obj for mapping ``state_id`` with ``state string``. The last entry in
            the dictionary should be the "unknown" state.

        assign_file : str, default : None
            A string pointing to the ``assign.h5`` file. Needed as a parameter for all functions,
            but is ignored if it's an MD trajectory.

        Returns
        -------
        dictionary : dict
            A dictionary mapping each ``state_id`` (float/int) with a ``state string`` (character).
            The last entry in the dictionary should be the "unknown" state.

        """
        # Other example for grouping multiple states into one.
        for idx, pathway in enumerate(data):
            # The following shows how you can "merge" multiple states into
            # a single one.
            pathway = numpy.asarray(pathway)
            # Further downsizing... to if pcoord is less than 5
            first_contact = numpy.where(pathway[:, 3] < 5)[0][0]
            for jdx, frame in enumerate(pathway):
                # First copy all columns over
                pathways[idx, jdx] = frame
                # ortho is assigned to state 0
                if frame[2] in [1, 3, 4, 6, 7, 9]:
                    frame[2] = 0
                # para is assigned to state 1
                elif frame[2] in [2, 5, 8]:
                    frame[2] = 1
                # Unknown state is assigned 2
                if jdx < first_contact:
                    frame[2] = 2
                pathways[idx, jdx] = frame

        # Generating a dictionary mapping each state
        dictionary = {0: 'A', 1: 'B', 2: '!'}

        return dictionary
