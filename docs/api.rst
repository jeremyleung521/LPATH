mPHAT API Documentation
=======================

.. autosummary::
   # :toctree: autosummary
   :recursive:

   mphat.discretize
   mphat.extract
   mphat.match
   mphat.mphat


Discretization
--------------
The discretize step allows you to assign MD trajectories (or WE simulations) into discrete state.

.. argparse::
   :module: mphat.argparser
   :func: add_discretize_args
   :prog: mphat

.. automodule:: mphat.discretize
    :members:


Extract
-------
The extract step allows you to extract successful trajectories from MD trajectories (or WE simulations).

.. automodule:: mphat.extract
    :members:

.. argparse::
   :module: mphat.argparser
   :func: add_extract_args
   :prog: mphat


.. automodule:: mphat.extract
    :only-members:

Match
-----
The match step allows you to compare and cluster pathways from the extract step.

.. argparse::
   :module: mphat.argparser
   :func: add_match_args
   :prog: mphat


.. automodule:: mphat.match
    :members:


mPHAT
-----
.. automodule:: mphat.mphat
    :members:


argparser
---------
.. automodule:: mphat.argparser
    :members:
