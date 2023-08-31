#!/bin/bash

# assign source and target states using the phi/psi angles from cMD simulation
lpath discretize --input "phi_psi.npy" --output "source_target.npy" --assign-function assign_states.assign_states

# extract all successful pathways connecting source and target states
# also extract the cluster labels for matching in the next step
lpath extract --extract-input "source_target.npy" --source-state 0 --target-state 1 --stride 1 --pcoord --extract-pcoord "states.npy" --feature-stride 1

# perform matching with condensing repeat pairs
# uses the cluster labels as states through reassign_custom
lpath match -ra reassign_custom.reassign_custom -op succ_traj/reassigned.pickle --condense 2
