#!/bin/bash

# assign source and target states using the phi/psi angles from WE simulations
lpath discretize -we -W ./multi.h5 --assign-arguments="--config-from-file --scheme C7_EQ -W multi.h5"

# extract all successful pathways connecting source and target states
# also extract the cluster labels for matching in the next step
lpath extract -we -W ./multi.h5 -A ./ANALYSIS/C7_EQ/assign.h5 -ss 0 -ts 1 -p -a labels --stride 25

# perform matching with condensing repeat pairs
# uses the cluster labels as states through reassign_custom
lpath match -we -ra reassign_custom.reassign_custom -op succ_traj/reassigned.pickle --condense 2
