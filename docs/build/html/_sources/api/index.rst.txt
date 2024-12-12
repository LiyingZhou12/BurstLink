.. BurstLink documentation master file, created by
   sphinx-quickstart on Mon Jun 24 03:12:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. module:: burstlink
.. automodule:: burstlink
   :noindex:

API Documentation
=================

Preprocessing: pp
-----------------

.. currentmodule:: burstlink.pp

.. autosummary::
   :toctree: .

   select_genepair_grn
   integration_grn
   RNAseq_analysis
   selection_GRNandRNAseq
   comparison


Tools: tl
---------

.. currentmodule:: burstlink.tl

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
   dge_analysis


Plotting: pl
------------

.. currentmodule:: burstlink.pl

.. autosummary::
   :toctree: .

   network_visualization
   TFs_interactiontype_network
   scatterplot_burst