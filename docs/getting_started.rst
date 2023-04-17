Getting Started
===============

This page details how to get started with ``mPHAT``.

1. Clone the repository from GitHub::

    git clone https://github.com/jeremyleung521/mphat


2. In the cloned directory, install the program with ``pip``. Specify ``westpa`` if you are analyzing a weighted ensemble simulation. Specify ``dev`` if you need to install all dependencies::

    cd mphat
    pip install -e .  # for essentials OR
    pip install -e .[westpa]  # with WESTPA OR
    pip install -e .[dev]  # all dependencies

3. Run all or parts of the program. See `Usage`_ for more information.


.. _Usage: https://mphat.readthedocs.io/en/latest/usage.html