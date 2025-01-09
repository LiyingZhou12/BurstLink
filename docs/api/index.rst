.. BurstLink documentation master file, created by
   sphinx-quickstart on Mon Jun 24 03:12:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


API Documentation
=================

Preprocessing
-------------

.. module:: burstlink.preprocessing
.. currentmodule:: burstlink

.. autosummary::
   :toctree: .

   select_genepair_grn
   integration_grn
   RNAseq_analysis
   selection_GRNandRNAseq
   comparison


Tools
-----

.. module:: burstlink.tools
.. currentmodule:: burstlink

.. autosummary::
   :toctree: .

   global_burst_link
   global_uni_burst_link
   genepair_inference
   genepair_burstinference
   genepair_interactionsinference
   ks_2samp
   tf_tg_analysis
   interaction_burst_regression
   burst_interaction_overall
   affinity_burst
   burst_interactionlevel_positive
   go_enrichment_analysis


Plotting
--------

.. module:: burstlink.plotting
.. currentmodule:: burstlink

.. autosummary::
   :toctree: .

   network_visualization
   TFs_interactiontype_network
   scatterplot_burst


Additional Resources
--------------------

.. toctree::
   :maxdepth: 2

   burstlink.preprocessing
   burstlink.tools
   burstlink.plotting
