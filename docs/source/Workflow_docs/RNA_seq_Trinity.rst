RNA-Seq without a reference genome, using Trinity 
-------------------------------------------------

:Author: Menachem Sklarz
:Affiliation: Bioinformatics Core Facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A pipeline for RNA-seq analysis using Trinity.
    
This workflow takes reads in `fastq` format, either paired-end or single, and assembles a trascriptome with `Trinity`.

It then runs ``align_and_estimate_abundance.pl`` and ``abundance_estimates_to_matrix.pl`` to map the reads to the trascriptome and create normalized counts tables. These tables can then be used in `DEseq2` or other tools for statistical analysis of RNA-seq data.
    
Steps:
~~~~~~~

1. Concatenating the read files into single files per direction (``merge``)

.. attention:: No QC steps are included here, but you should consider adding `trimmomatic` and `FastQC` steps to the workflow.

2. Adding tags required by trinity to the read titles (/1 and /2 for F and R. See `Running-Trinity <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity>`_.
3. Running ``Trinity``. Trinity must be configured to run on a cluster. The configuration file is passed with ``--grid_conf``.
4. Mapping of the reads is performed with ``trinity_mapping`` module.
5. Creating statistical tables is performed with ``trinmap_statistics`` module.
    
Workflow Schema
~~~~~~~~~~~~~~~~

.. image:: RNA_seq_Trinity.png
   :alt: Denovo RNA seq DAG

Requires
~~~~~~~~

`fastq` files. Paired end or single-end.

Programs required
~~~~~~~~~~~~~~~~~~

* `bowtie2      <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
* `Trinity      <https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_
* `RSEM         <https://deweylab.github.io/RSEM/>`_
* `samtools     <http://www.htslib.org/>`_
* `STAR         <https://github.com/alexdobin/STAR>`_


Example of Sample File
~~~~~~~~~~~~~~~~~~~~~~

::

    Title	RNA_seq_denovo

    #SampleID	Type	Path    lane
    Sample1	Forward	/path/to/Sample1_F1.fastq.gz 1
    Sample1	Forward	/path/to/Sample1_F2.fastq.gz 2
    Sample1	Reverse	/path/to/Sample1_R1.fastq.gz 1
    Sample1	Reverse	/path/to/Sample1_R2.fastq.gz 2
    Sample2	Forward	/path/to/Sample2_F1.fastq.gz 1
    Sample2	Reverse	/path/to/Sample2_R1.fastq.gz 1
    Sample2	Forward	/path/to/Sample2_F2.fastq.gz 2
    Sample2	Reverse	/path/to/Sample2_R2.fastq.gz 2

Download
~~~~~~~~~

The workflow file is available :download:`here <../../../Workflows/RNA_seq_Trinity.yaml>`

