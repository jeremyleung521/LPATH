![lpath_logo.png](lpath_logo.png)

Linguistics Pathway Analysis of Trajectories using Hierarchical clustering
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/jeremyleung521/mphat/workflows/CI/badge.svg)](https://github.com/jeremyleung521/mphat/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/mphat/badge/?version=latest)](https://mphat.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/jeremyleung521/mphat/branch/main/graph/badge.svg?token=5SS08RH1MO)](https://codecov.io/gh/jeremyleung521/mphat)

A user-friendly, Python tool for clustering pathways from molecular dynamics and weighted ensemble simulations.

### Copyright

Copyright (c) 2023, Anthony Bogetti, Jeremy Leung, Lillian Chong

### TODO
- [x] EXTRACT step only does WE right now, need to extend to plain MD (i.e. a states.npy)
- [x] Usage.rst page
- [ ] example ipynbs are bare bones, need sample files
- [ ] unit testing
- [x] extract defaults to last frame of iteration. Implement stride and actually give pcoord/auxdata to correct frame.
- [x] command line interface/argparser
- [x] Sphinx/Read the Docs autodoc
- [x] Run each tool directly
- [ ] Rename Repo

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
