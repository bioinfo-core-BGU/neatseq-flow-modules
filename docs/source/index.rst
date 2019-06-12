.. neatseq_flow modules documentation master file, created by
   sphinx-quickstart on Sun Jan 08 15:32:48 2017.

.. _module_repo_docs:

**NeatSeq-Flow**'s Module and Workflow Repository 
===================================================

.. .. figure:: figs/NeatSeq_Flow_logo.png
   :scale: 60 %
   :align: center
   :alt: NeatSeq-Flow logo
   :target: https://neatseq-flow.readthedocs.io/en/latest/index.html#

.. image:: https://readthedocs.org/projects/neatseq-flow-modules/badge/?version=latest
   :target: http://neatseq-flow-modules.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

This repository contains modules and workflows for use with **NeatSeq-Flow**.  To install, follow the instructions in :ref:`installation section bellow <installation>`.

See also the `main NeatSeq-Flow website <https://neatseq-flow.readthedocs.io/en/latest/>`_.



.. _neatseq_flow_workflows:

.. toctree::
   :maxdepth: 1
   :caption: Workflows

   Workflow_docs/Tutorial
   Workflow_docs/RNA_seq_reference
   Workflow_docs/RNA_seq_Trinity
   Workflow_docs/QIIME_workflow
   Workflow_docs/QIIME2_workflow
   Workflow_docs/GATK_workflow
   Workflow_docs/Microbe-Flow
   Workflow_docs/ChIP_seq
   Workflow_docs/Metagenomics

..   Workflow_docs/Clustering
..   Workflow_docs/BLAST_db
..   Workflow_docs/BLAST_fasta
..   Workflow_docs/Assembly_Indexing_mapping


.. _neatseq_flow_modules:

.. toctree::
   :maxdepth: 2
   :caption: Modules

   Module_docs/AllModules
..   Module_docs/PreparationAndQC
..   Module_docs/Mapping
..   Module_docs/BAMConversion
..   Module_docs/ChIPseq
..   Module_docs/GenomeAssembly
..   Module_docs/TranscriptomeAssembly
..   Module_docs/TranscriptomeAnnotation
..   Module_docs/SequenceAnnotation
..   Module_docs/SequenceSearching
..   Module_docs/Variants
..   Module_docs/Metagenomics
..   Module_docs/Microbiology
..   Module_docs/QIIME
..   Module_docs/QIIME2
..   Module_docs/GATK
..   Module_docs/SequenceClustering
..   Module_docs/VariousReportingPrograms
..   Module_docs/Miscellaneous
..   Module_docs/GenericModules


..   all_modules
   

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   Module_docs/Generic_module






.. _installation:

Installing and using the modules
-----------------------------------

**NeatSeq-Flow** is available for download on `github <https://github.com/bioinfo-core-BGU/neatseq-flow>`_.


.. attention::
   In order to include these modules in your workflow, please:

   #. download the repository::

       curl -LO https://github.com/bioinfo-core-BGU/neatseq-flow-modules/archive/master.zip

   #. **Alternatively**, clone the repository::

       git clone https://github.com/bioinfo-core-BGU/neatseq-flow-modules.git

   #. add the following line in the `Global_params` section of your workflow parameter file::

       module_path: /path/to/location/of/modules/repository

   See more about definition of workflow parameter files in the NeatSeq-Flow User Manual at the `Parameter file definition <http://neatseq-flow.readthedocs.io/en/latest/02a.FileDefinition.html#parameter-file-definition>`_ section.

.. Note::
   Some of the modules in this package are included in the main **NeatSeq-Flow** repository. These are indicated below with a :sup:`*`.

.. important::
   **NeatSeq-Flow** enables users to program their own modules and workflows. You are encouraged to share your modules with the public by adding it to this repository. In order to do so, please fork the repository on github, upload your new module or workflow and open a pull request.
