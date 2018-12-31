Microbiome analysis using QIIME2
--------------------------------

:Author: Menachem Sklarz
:Affiliation: Bioinformatics Core Facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A workflow for executing the `Moving Windows tutorial <https://docs.qiime2.org/2018.11/tutorials/moving-pictures/#moving-pictures-tutorial>`_ with QIIME2.


Steps:
~~~~~~~

#. merge
#. import sequence data
#. Demultiplex and show statistics
#. `dada2 <https://benjjneb.github.io/dada2/>`_  and visualization.
#. Building a phylogenetic tree
#. Core diversity analysis
#. Comparing alpha and beta groups differences.
#. Creating emperor visualizations.
#. Creating &alpha;-rarefaction curves.
#. Taxonomic classification and visualizations.


Workflow Schema
~~~~~~~~~~~~~~~~

.. image:: QIIME2_workflow.jpg
   :alt: QIIME2 workflow DAG

Requires
~~~~~~~~

* `Sample metadata <https://docs.qiime2.org/2018.11/tutorials/moving-pictures/#sample-metadata>`_
* barcodes file
* sequence file

The latter two should be downloaded and stored in a directory called ``emp-single-end-sequences`` as described in the `Obtaining and importing data <https://docs.qiime2.org/2018.11/tutorials/moving-pictures/#obtaining-and-importing-data>`_ section of the `tutorial webpage <https://docs.qiime2.org/2018.11/tutorials/moving-pictures/>`_.

Programs required
~~~~~~~~~~~~~~~~~~

* `QIIME2       <https://qiime2.org/>`_


Example of Sample File
~~~~~~~~~~~~~~~~~~~~~~

::

   Title	MovingPictures


   #Type	Path
   EMPSingleEndSequences	/path/to/emp-single-end-sequences
   metadata	/path/to/sample-metadata.tsv
   TaxonomicClassifier	https://data.qiime2.org/2018.11/common/gg-13-8-99-515-806-nb-classifier.qza

Download
~~~~~~~~~

The workflow file is available :download:`here <../../../Workflows/QIIME2_MovingPic.yaml>`

