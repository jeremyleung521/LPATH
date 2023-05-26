Getting Started
===============

This page details how to get started with ``lpath``.

1. Clone the repository from GitHub::

    git clone https://github.com/jeremyleung521/lpath

2. In the cloned directory, install the program with ``pip``. Specify ``westpa`` if you are analyzing a weighted ensemble simulation. Specify ``dev`` if you need to install all dependencies::

    cd lpath
    pip install -e .  # for essentials OR
    pip install -e .[westpa]  # with WESTPA OR
    pip install -e .[dev]  # all dependencies

3. Run all or parts of the program. See `Usage`_ for more information.


.. _Usage: https://lpath.readthedocs.io/en/latest/usage.html