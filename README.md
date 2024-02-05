![lpath_logo.png](https://raw.githubusercontent.com/chonglab-pitt/LPATH/main/logo/lpath_logo.png)

Linguistics Pathway Analysis of Trajectories using Hierarchical clustering
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/chonglab-pitt/lpath/workflows/CI/badge.svg)](https://github.com/chonglab-pitt/lpath/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/lpath/badge/?version=latest)](https://lpath.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/624926339.svg)](https://zenodo.org/badge/latestdoi/624926339)
[![Publication](https://img.shields.io/badge/Publication-darkred)](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c01318)

A user-friendly, Python tool for clustering pathways from molecular dynamics and weighted ensemble simulations.


## Copyright

This software is distributed with the MIT License.

Copyright © 2023, Anthony Bogetti, Jeremy Leung, Lillian Chong


## Quickstart Guide

Users may install this program via [PyPI](https://pypi.org/project/lpath-md/) (recommended) or 
from [GitHub Source](https://github.com/chonglab-pitt/lpath). By default, LPATH is installed with the bare minimum to 
run analysis on classic MD simulations. Users who intend to analyze WESTPA simulations should install WESTPA 
in the same environment.

1. Install from [PyPI](https://pypi.org/project/lpath-md/):
    ```python -m  pip install lpath-md```
2. Install from [GitHub Source](https://github.com/chonglab-pitt/lpath):
    1. Clone the git repository:  ```git clone https://github.com:chonglab-pitt/lpath```
    2. Move into the folder: ```cd lpath```
    3. Install the program in editable mode: ```python -m pip install -e .```
3. There are a variety of optional dependencies you may install. These work with either installation options listed above.
    1. Install with `WESTPA`: ```python -m pip install lpath-md[westpa]```
    2. Install with terminal user interface (TUI): ```python -m pip install lpath-md[tui]```
    3. Developers can install with [dev] (all dependencies) or [test] (minimal dependencies to run tests): ```python -m pip install lpath-md[dev]```
    4. Options may be combined: ```python -m pip install lpath-md[westpa,tui]```


### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
