Common Errors
=============

This page details common errors you might encounter while using ``lpath``.

Discretize
__________

[NONE AT THE MOMENT]

Extract
_______

[NONE AT THE MOMENT]

Match
_____

IndexError and DivideByZeroError

The following tracebacks imply that your states have all accidentally been cleared out by the ``expand_shorter_traj()`` function::
    Traceback(most recent call last):
      File "/home/user/apps/miniconda3/envs/westpa-lpath/bin/lpath", line 8, in <module>
        sys.exit(entry_point())
                 ^^^^^^^^^^^^^
      File "/home/user/Documents/lpath/lpath/lpath.py", line 52, in entry_point
        args.func(args)
      File "/home/user/Documents/lpath/lpath/match.py", line 492, in main
        dist_matrix, weights = gen_dist_matrix(pathways, dictionary, file_name=arguments.dmatrix_save,
                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/Documents/lpath/lpath/match.py", line 280, in gen_dist_matrix
        weights.append(nonzero[-1][-1])
                       ~~~~~~~^^^^
      IndexError: index -1 is out of bounds for axis 0 with size 0


    Traceback (most recent call last):
      File "/home/atb43/apps/miniconda3/envs/westpa2-ibstates-fix/bin/lpath", line 8, in <module>
        sys.exit(entry_point())
                 ^^^^^^^^^^^^^
      File "/home/atb43/Documents/lpath-stride/lpath/lpath.py", line 52, in entry_point
        args.func(args)
      File "/home/atb43/Documents/lpath-stride/lpath/match.py", line 601, in main
        dist_matrix, weights = gen_dist_matrix(pathways, dictionary, file_name=arguments.dmatrix_save,
                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/Documents/lpath-stride/lpath/match.py", line 392, in gen_dist_matrix
        distmat = pairwise_distances(
                  ^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/apps/miniconda3/envs/westpa2-ibstates-fix/lib/python3.11/site-packages/sklearn/metrics/pairwise.py", line 2039, in pairwise_distances
        return _parallel_pairwise(X, Y, func, n_jobs, **kwds)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/apps/miniconda3/envs/westpa2-ibstates-fix/lib/python3.11/site-packages/sklearn/metrics/pairwise.py", line 1579, in _parallel_pairwise
        return func(X, Y, **kwds)
               ^^^^^^^^^^^^^^^^^^
      File "/home/atb43/apps/miniconda3/envs/westpa2-ibstates-fix/lib/python3.11/site-packages/sklearn/metrics/pairwise.py", line 1623, in _pairwise_callable
        out[i, j] = metric(X[i], Y[j], **kwds)
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/Documents/lpath-stride/lpath/match.py", line 393, in <lambda>
        X=path_strings, metric=lambda x, y: calc_dist(x, y, dictionary)
                                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      File "/home/atb43/Documents/lpath-stride/lpath/match.py", line 69, in calc_dist
        similarity = (2 * km) / (int(len(seq1_str) + len(seq2_str)))
                     ~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ZeroDivisionError: division by zero
    

There are two common cause for this problem due to the custom reassign function:
    1. The most common cause is that you have utilized a custom reassign function but forgot to include an "unknown state" in the dictionary. By LPATH conventions, the last entry in the dictionary is always reserved for this state, represented by ``!`` (or ``chr(33)``). During string comparison, any states assigned "unknown" are removed. Adding the following line into your custom reassign function will typically solve it::
        dictionary[len(dictionary)] = '!'

    2. Another potential cause is that you forgot to copy the other columns of data (iteration, seg_id etc.) into the ``pathways`` array during reassignment. The ``expand_shorter_traj()`` function assumes all segments with iteration 0 are there for padding (i.e., there because the length of the transition is < the longest successful trajectory extracted) and are assigned to the "unknown state". Include the following line within the second ``for`` loop to copy all columns before reassignment::
        for idx, val2 in enumerate(flipped_val):
            pathways[idx,idx2] = val2

Plot
----

[NONE AT THE MOMENT]