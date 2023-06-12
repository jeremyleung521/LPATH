LPATH API Documentation
=======================

.. autosummary::
   # :toctree: autosummary
   :recursive:

   lpath.discretize
   lpath.extract
   lpath.match
   lpath.plot
   lpath.lpath


Common Parameters
-----------------
All steps take the following parameters.

.. argparse::
   :module: lpath.argparser
   :func: add_common_args
   :prog: lpath {discretize,extract,match,plot,all}


Discretization Step
-------------------
The ``discretize`` step allows you to assign MD trajectories (or WE simulations) into discrete states.

.. argparse::
   :module: lpath.argparser
   :func: add_discretize_args
   :prog: lpath discretize

.. automodule:: lpath.discretize
    :members:


Extract Step
------------
The ``extract`` step allows you to extract successful trajectories from MD trajectories (or WE simulations) based on definitions from ``disretize``.

.. argparse::
   :module: lpath.argparser
   :func: add_extract_args
   :prog: lpath extract

.. automodule:: lpath.extract
    :members:


Match Step
----------
The ``match`` step allows you to compare and cluster pathways from the ``extract`` step.

.. argparse::
   :module: lpath.argparser
   :func: add_match_args
   :prog: lpath match


.. automodule:: lpath.match
    :members:


Match Step
----------
The ``plot`` step allows you to plot clusters and pathways from the ``match`` step.

.. argparse::
   :module: lpath.argparser
   :func: add_plot_args
   :prog: lpath plot


.. automodule:: lpath.plot
    :members:


All LPATH Steps
---------------
.. argparse::
   :module: lpath.argparser
   :func: add_all_args
   :prog: lpath all

.. automodule:: lpath.main
    :members:


argparser
---------
.. automodule:: lpath.argparser
    :members:
