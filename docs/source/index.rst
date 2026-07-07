MOT-MIC: multi-organ optimal transport for metastasis-initiating cell discovery
================================================================================

.. image:: _static/img/motmic_head_image.png
   :width: 520px
   :align: left

MOT-MIC is a computational framework for discovering metastasis-initiating
cells from paired primary and metastatic single-cell transcriptomes.

It provides a unified workflow for scoring primary tumor cells by metastatic
competence, assigning organ-specific metastatic propensity, validating
predictions against lineage and clinical evidence, and ranking metastatic genes
with SHAP-based interpretability.

.. raw:: html

   <br clear="left"/>

Key features
------------

- Multi-organ MIC scoring from paired primary and metastatic scRNA-seq data.
- Unbalanced optimal transport with sparse top-k origin filtering.
- Organotropic MIC scores for liver, lung, bone, brain, or user-defined sites.
- Lineage-aware validation with GSE173958.
- scTour-style tutorials for discovery, lineage validation, clinical validation,
  and spatial transfer.
- SHAP-based prioritization of pan-MIC and organ-specific metastatic genes.

Installation
------------

MOT-MIC requires Python 3.8 or later::

    git clone https://github.com/yuhan1li/scMIC.git
    cd scMIC
    pip install -r requirements.txt

Run the quick example::

    python scripts/run_example.py

Tutorials
---------

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   :hidden:

   notebook/01_motmic_example
   notebook/02_GSE173958_lineage_validation
   notebook/03_GSE249057_timecourse_discovery
   notebook/04_GSE178318_human_crc_liver_metastasis
   notebook/05_GSE277783_spatial_validation

API
---

.. toctree::
   :maxdepth: 2
   :caption: API
   :hidden:

   api_core
   api_io
   api_interpret
   api_validation

Reference datasets
------------------

The recommended validation hierarchy is GSE173958 lineage tracing, GSE249057
time-course discovery, GSE178318 paired human CRC liver metastasis, GSE277783
spatial validation, and TCGA bulk survival transfer.

