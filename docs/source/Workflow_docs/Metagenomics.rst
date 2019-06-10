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
    1. Assembly is done per-sample with ``spades``.
    2. The assemblies are quality-tested with ``quast``.
    3. Assemblies are annotated with ``Prokka``.
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
* `prokka       <http://www.vicbioinformatics.com/software.prokka.shtml>`_
* `quast        <http://bioinf.spbau.ru/quast>`_
* `spades       <http://bioinf.spbau.ru/spades>`_
* `RGI          <https://card.mcmaster.ca/analyze/rgi>`_
* `KronaTools   <https://github.com/marbl/Krona/wiki/KronaTools>`_

All the programs used in this workflow can be installed with conda. See section :ref:`quick-conda-start`_ below.

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


.. _quick-conda-start:

Quick start with conda
~~~~~~~~~~~~~~~~~~~~~~~

For easy setup of the workflow, including a sample dataset, use the following instructions for complete installation with conda:

#. Create and activate a conda environment with all the required programs::

    wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/docs/source/Workflow_docs/Metagenomics_conda.yaml
    conda env create -f Metagenomics_conda.yaml
    conda activate Metagenomics

#. Create a sample file. It should look like the following, only the file names should be replaced with absolute file names::

        Title   Trinity_example

        #SampleID       Type    Path
        Sample1 Forward 00.Raw_reads/reads.left.fq.gz
        Sample1 Reverse 00.Raw_reads/reads.right.fq.gz

   .. Tip:: To get the full path to a file, use the following command:

      .. code-block:: bash

         readlink -f 00.Raw_reads/reads.left.fq.gz

#. Create a directory for your databases. Save the location of the directory in $DBDIR.

   .. code-block:: bash

     export DBDIR=/path/to/databases_dir
     mkdir -p $DBDIR

#. **Install required databases**

    .. Tip:: File ``Metagenomics_DBinstall_cmds.sh`` contains a script for installing all the databases described below.

       Execution might take a while due to the large datasetb being downloaded, therefore it is recommended to execute as follows (**After setting $DBDIR!!!**):

       .. code-block:: bash

          wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/docs/source/Workflow_docs/Metagenomics_DBinstall_cmds.sh
          nohup bash Metagenomics_DBinstall_cmds.sh &

    #. metaphlan:

       Running metaphlan will download the database for you:

       .. code-block:: bash

            metaphlan2.py \
                --input_type fastq \
                --bowtie2_exe bowtie2 \
                --bowtie2db $DBDIR/MetaPhlAn_temp


    #. kraken2:

       Installing kraken2 database takes a long time and requires a considerable amount of disk space.

       .. code-block:: bash

            mkdir -p $DBDIR/kraken2
            kraken2-build \
                --standard \
                --threads 10 \
                --db $DBDIR/kraken2

       .. Attention::  If ``rsync`` dosen't work for you, you can try adding the ``--use-ftp`` to the ``kraken2-build`` command to use ``wget`` instead.

    #. centrifuge:

       .. code-block:: bash

            mkdir -p $DBDIR/centrifuge
            centrifuge-download \
                -o $DBDIR/centrifuge/taxonomy \
                taxonomy

            centrifuge-download \
                -o $DBDIR/centrifuge \
                -m -d "archaea,bacteria,viral" refseq \
                > $DBDIR/centrifuge/seqid2taxid.map

            cat $DBDIR/centrifuge/*/*.fna > $DBDIR/centrifuge/input-sequences.fna

            mkdir $DBDIR/centrifuge/index
            centrifuge-build -p 4 \
                --conversion-table $DBDIR/centrifuge/seqid2taxid.map \
                --taxonomy-tree $DBDIR/centrifuge/taxonomy/nodes.dmp \
                --name-table $DBDIR/centrifuge/taxonomy/names.dmp \
                $DBDIR/centrifuge/input-sequences.fna \
                $DBDIR/centrifuge/index/arch_bac_vir


        .. Attention:: The download commands may fail because of the libssl version.

    #. krona:

       .. code-block:: bash

            ktUpdateTaxonomy.sh $DBDIR/krona/taxonomy

    #. kaiju:

       Kaiju provides different databases which can be downloaded. To get a list of options, just execute ``kaiju-makedb`` with no arguments:

       The following commands demonstrate how to get the ``nr`` database including eukaryotes (``nr_euk``) and the ``progenomes`` database.

       .. code-block:: bash

            mkdir -p $DBDIR/kaiju
            cd $DBDIR/kaiju
            kaiju-makedb -s progenomes -t 10
            kaiju-makedb -s nr_euk -t 10
            cd -

    #. HUMAnN2:

       `Online help on downloading databases <https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-5-download-the-databases>`_.

       .. code-block:: bash

            mkdir -p databases/HUMAnN2
            humann2_databases --download chocophlan full  $DBDIR/HUMAnN2
            humann2_databases --download uniref uniref90_diamond  $DBDIR/HUMAnN2/uniref90
            humann2_databases --download uniref uniref50_diamond  $DBDIR/HUMAnN2/uniref50

            humann2_config --update database_folders nucleotide $DBDIR/HUMAnN2/chocophlan
            humann2_config --update database_folders protein $DBDIR/HUMAnN2/uniref90

       .. Attention:: The commands download the recommended translated databases. For other options, see
            the `Download a translated search database <https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-download-a-translated-search-database>`_ section of the tutorial.

#. Get the parameter file with::

    wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/Workflows/Menagenomics.yaml

#. **Settings to set in the parameter file**

   You will have to make some changes to the parameter file to suit your needs:

   #. Set the parameters in the ``Global_params`` section to suit your cluster. Alternatively, set ``Executor`` to ``Local`` for running on a single machine.
   #. In the ``Vars`` section, set ``database_prefix`` to the location of your databases dir, which is the value of ``$DBDIR`` set above. If $DBDIR is set, you can use the following ``sed`` command to set the ``database_prefix`` correctly:

      .. code-block:: bash

         sed -i s+\$DBdir+$DBDIR+ Metagenomics12.yaml

   #. In ``Vars.databases.kaiju``, you will have to make sure the value of ``fmi`` fits the database you decide to use. In the provided parameter file, the ``nr_euk`` is set. The equivalent ``fmi`` value for the ``progenomes`` database is commented out.
   #. Go over the ``redirects`` sections in the parameter file and make sure they are set according to your requirements.
   #. If you have a fasta file with sequences to search for within your metagenome assemblies, set the ``proteins_of_interest`` variable to the full path to that file. If not, you can delete or uncomment the ``SKIP`` line in steps ``make_blast_db_per_assembly``, ``blast_proteins_vs_assemblies`` and ``parse_blast``.


#. In the conda definitions (line 46), set ``base:`` to the path to the conda installation which you used to install the environment.

    You can get the path by executing the following command, **when inside the Metagenomics conda environment**:

    .. code-block:: bash

        echo $CONDA_EXE | sed -e 's/\/bin\/conda$//g'

#. `Execute NeatSeq-Flow  <https://neatseq-flow.readthedocs.io/en/latest/02b.execution.html#executing-neatseq-flow>`_.


.. Tip:: See also `this nice presentation <http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2015/07/cfar_lab_09182015.pdf>`_ by Galeb Abu-Ali, Eric Franzosa and Curtis Huttenhower


