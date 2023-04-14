mPHAT API Documentation
=======================

.. autosummary::
   # :toctree: autosummary
   :recursive:

   mphat.discretize
   mphat.extract
   mphat.match
   mphat.mphat


Discretization Step
-------------------
The ``discretize`` step allows you to assign MD trajectories (or WE simulations) into discrete states.

.. argparse::
   :module: mphat.argparser
   :func: add_discretize_args
   :prog: mphat

.. automodule:: mphat.discretize
    :members:


Extract Step
------------
The ``extract`` step allows you to extract successful trajectories from MD trajectories (or WE simulations) based on definitions from ``disretize``.

.. argparse::
   :module: mphat.argparser
   :func: add_extract_args
   :prog: mphat


.. automodule:: mphat.extract
    :members:


Match Step
----------
The ``match`` step allows you to compare and cluster pathways from the ``extract`` step.

.. argparse::
   :module: mphat.argparser
   :func: add_match_args
   :prog: mphat


.. automodule:: mphat.match
    :members:


All mPHAT Steps
---------------
.. automodule:: mphat.mphat
    :members:


argparser
---------
.. automodule:: mphat.argparser
    :members:
