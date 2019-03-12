Shotgun metagenomics using kraken, MetaPhlAn2, Kaiju and HUMAnN2
----------------------------------------------------------------

:Author: Menachem Sklarz
:Affiliation: Bioinformatics Core Facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A workflow for executing various analyses on metagenomics data.

The workflow uses two approaches:

1. Analysis of the raw reads.
2. Assembly of the reads and analysis of the assembled contigs

**Developed as part of a study led by Prof. Jacob Moran-Gilad.**
 
Steps:
~~~~~~~

a. Analysis of the raw reads with:
    * ``kraken``
    * ``metaphlan2``
    * ``kaiju``
    * ``HUMAnN2``

    The output from the former three programs is also plotted with krona (to disable plotting with krona, comment out the lines referring to krona in the instance definition.)  
b. Assembly and analysis of the assembled reads:
    1. Assembly is done with two tools: ``spades`` and ``megahit``.
    2. Each assembly is quality tested with ``quast``.
    3. assemblies are annotated with ``Prokka``.
    4. Antibiotic resistance is determined with ``CARD_RGI``.
    5. **Not included**. Resistance and virulence can also be determined by BLASTing AR and virulence databases against the assemblies. See module BLAST.

Workflow Schema
~~~~~~~~~~~~~~~~

.. image:: Metagenomics.png   
   :alt: Metagenomics DAG

Requires
~~~~~~~~

`fastq` files. Paired end or single-end.

Programs required
~~~~~~~~~~~~~~~~~~

* `FastQC       <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
* `trimmomatic  <http://www.usadellab.org/cms/?page=trimmomatic>`_
* `kraken       <https://ccb.jhu.edu/software/kraken/>`_
* `kaiju        <http://kaiju.binf.ku.dk/>`_
* `metaphlan2   <https://bitbucket.org/biobakery/metaphlan2>`_
* `bowtie2      <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
* `diamond      <https://ab.inf.uni-tuebingen.de/software/diamond>`_
* `HUMAnN2      <http://huttenhower.sph.harvard.edu/humann2>`_
* `megahit      <https://github.com/voutcn/megahit>`_
* `prokka       <http://www.vicbioinformatics.com/software.prokka.shtml>`_
* `quast        <http://bioinf.spbau.ru/quast>`_
* `spades       <http://bioinf.spbau.ru/spades>`_
* `RGI          <https://card.mcmaster.ca/analyze/rgi>`_
* `KronaTools   <https://github.com/marbl/Krona/wiki/KronaTools>`_


Example of Sample File
~~~~~~~~~~~~~~~~~~~~~~

::

    Title	Metagenomics

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

The workflow file is available :download:`here <../../../Workflows/Metagenomics.yaml>`


Quick start with conda
~~~~~~~~~~~~~~~~~~~~~~~

For easy setup of the workflow, including a sample dataset, use the following instructions for complete installation with conda:

.. Attention:: ``rnammer`` is not available with CONDA. To use it, you will have to install it and modify it `following the instructions here <https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required#rnammer-free-academic-download>`_.

#. Create and activate a conda environment with all the required programs::

    wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/docs/source/Workflow_docs/Metagenomics_conda.yaml
    conda env create -f RNA_seq_Trinity_conda.yaml
    conda activate RNA_trinity

#. Get the raw data from Trinity::

    mkdir 00.Raw_reads
    cp $CONDA_PREFIX/opt/trinity-2.8.4/Docker/test_data/reads.right.fq.gz 00.Raw_reads/
    cp $CONDA_PREFIX/opt/trinity-2.8.4/Docker/test_data/reads.left.fq.gz 00.Raw_reads/

#. Create a sample file. It should look like the following, only the file names should be replaced with absolute file names::

        Title   Trinity_example

        #SampleID       Type    Path
        Sample1 Forward 00.Raw_reads/reads.left.fq.gz
        Sample1 Reverse 00.Raw_reads/reads.right.fq.gz

   .. Tip:: To get the full path to a file, use the following command:

      .. code-block:: bash

         readlink -f 00.Raw_reads/reads.left.fq.gz

#. Get the parameter file with::

    wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/Workflows/RNA_seq_Trinity.yaml

#. In the conda definitions (line 46), set ``base:`` to the path to the conda installation which you used to install the environment.

    You can get the path by executing the following command::

        echo $CONDA_EXE | sed -e 's/\/bin\/conda$//g'


#. If you want to use Trinotate, create a directory for the required databases (this step takes some time to complete)::

    mkdir Trinotate_dbs;
    Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate_dbs/Trinotate

    mv uniprot_sprot.* Trinotate_dbs/
    mv Pfam-A.hmm.gz Trinotate_dbs/
    cd Trinotate_dbs/
    makeblastdb -in uniprot_sprot.pep -dbtype prot
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    cd -

.. Attention:: If you already have the Trinotate databases downloaded and setup, you do not have to do the last steps. You can set the paths to the databases in the ``databases`` subsection of the ``Vars`` section in the parameter file.

#. If you want to use BUSCO:

    #. Download a template config file with the following command and edit is as necessary::

        wget -O config.ini https://gitlab.com/ezlab/busco/raw/master/config/config.ini.default

    #. Set the Vars.databases.BUSCO variable to the URL or the BUSCO dataset to use. Choose a URL from this list: `<https://busco.ezlab.org/frame_wget.html>`_.

#. `Execute NeatSeq-Flow  <https://neatseq-flow.readthedocs.io/en/latest/02b.execution.html#executing-neatseq-flow>`_.
