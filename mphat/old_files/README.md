# ChongLabPrivate/analysis-scripts/pathway-analysis
A folder containing pathway analysis scripts written by Anthony Bogetti and Jeremy Leung (2022).

Basic workflow:
  1.  If working with multiple replicates of WE simulation, combine with `w_multi_west` (along with `--ibstates` flag for use with westpa.analysis)
  2.  Assign each frame to a state using `w_assign`.
  3.  `01_gen_succ_list.py` (and the equivalent ray version `01_gen_succ_list_ray.py`) goes through all the states and outputs all "successful trajectories" in a pickle object.
  4. `02_pattern_match_v6.py` will uses a string comparison algorithm to determine similarity of the sequence of states traversed. Dendrogram will be generated to help determine optimal number of clusters.
