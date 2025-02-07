Getting Started
===============

This page details how to get started with ``lpath``.

Users may install this program via `PyPI`_ (recommended) or
from `GitHub Source`_. By default, LPATH is installed with the bare minimum to
run analysis on classic MD simulations. Users who intend to analyze WESTPA simulations should install WESTPA 
in the same environment.

1. Install from `PyPI`_:

    python -m  pip install lpath

2. Install from `GitHub Source`_ :
    1. Clone the git repository:  ``git clone https://github.com:chonglab-pitt/lpath``
    2. Move into the folder: ```cd lpath```
    3. Install the program in editable mode: ``python -m pip install -e .``
3. There are a variety of optional dependencies you may install. These work with either installation options listed above.
    1. Install with ``WESTPA``: ``python -m pip install lpath[westpa]``
    2. Install with terminal user interface (TUI): ``python -m pip install lpath[tui]``
    3. Developers can install with [dev] (all dependencies) or [test] (minimal dependencies to run tests): ```python -m pip install lpath[dev]```
    4. Options may be combined: ``python -m pip install lpath[westpa,tui]``
4. Run all or parts of the program. See `Usage`_ for more information.

.. _PyPI: https://pypi.org/project/lpath/
.. _GitHub Source: https://github.com/chonglab-pitt/lpath
.. _Usage: https://lpath.readthedocs.io/en/latest/usage.html
