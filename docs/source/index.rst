.. neatseq_flow modules documentation master file, created by
   sphinx-quickstart on Sun Jan 08 15:32:48 2017.

.. _module_repo_docs:

NeatSeq_Flow's Module and Workflow Repository 
=============================================

.. figure:: figs/NeatSeq_Flow_logo.png
   :scale: 60 %
   :align: center
   :alt: NeatSeq-Flow logo

.. image:: https://readthedocs.org/projects/neatseq-flow-modules/badge/?version=latest
   :target: http://neatseq-flow-modules.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

                
**NeatSeq-Flow** is a lightweight software for efficient execution of high throughput sequencing workflows. See the **NeatSeq-Flow** documentation at http://neatseq-flow.readthedocs.io/en/latest/.

**NeatSeq-Flow** is available for download on `github <https://github.com/bioinfo-core-BGU/neatseq-flow>`_.

Following is a list of modules included in this repository.

.. attention:: 
   In order to include these modules in your workflow, please:
   
   1. download the repository::
   
       wget https://github.com/bioinfo-core-BGU/neatseq-flow-modules/archive/master.zip

   2. **Alternatively**, clone the repository::
   
       git clone https://github.com/bioinfo-core-BGU/neatseq-flow-modules.git
      
   3. add the following line in the `Global_params` section of your workflow parameter file::
   
       module_path: /path/to/location/of/modules/repository
   
   See more about definition of workflow parameter files in the documentation for the main **NeatSeq-Flow** package at http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#parameter-file-definition.
   
.. Note::
   Some of the modules in this package are included in the main **NeatSeq-Flow** repository. These are indicated below with a :sup:`*`.
   
.. important::
   **NeatSeq-Flow** enables users to program their own modules and workflows. You are encouraged to share your modules with the public by adding it to this repository. In order to do so, please fork the repository on github, upload your new module or workflow and open a pull request.


.. Modules included in the repository 
.. -----------------------------------


.. toctree::
   :maxdepth: 2
   :caption: Modules

   Module_docs/PreparationAndQC
   Module_docs/Mapping
   Module_docs/BAMConversion
   Module_docs/ChIPseq
   Module_docs/GenomeAssembly
   Module_docs/TranscriptomeAssembly
   Module_docs/SequenceAnnotation
   Module_docs/SequenceSearching
   Module_docs/Variants
   Module_docs/Metagenomics
   Module_docs/Microbiology
   Module_docs/QIIME
   Module_docs/SequenceClustering
   Module_docs/VariousReportingPrograms
   Module_docs/GenericModule

..   all_modules
   
   
   
.. Workflows included in the repository
.. -------------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Workflows

   Workflow_docs/Assembly_Indexing_mapping
   Workflow_docs/BLAST_db
   Workflow_docs/BLAST_fasta
   Workflow_docs/ChIP_seq_bowtie2
   Workflow_docs/Metagenomics
   Workflow_docs/RNA_seq_Trinity
   Workflow_docs/RNA_seq_reference
   Workflow_docs/variant_calling
   Workflow_docs/Clustering
   Workflow_docs/QIIME_workflow


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

