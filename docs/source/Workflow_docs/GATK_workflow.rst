GATK variant analysis
----------------------

A workflow for somatic short variant discovery (SNVs + Indels) (based on `GATK Best Practices <https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146>`_. )

 
Steps:
~~~~~~~

1. Read preparation:
    1. merge
    2. trimmomatic - For cleaning the reads
    3. FastQC - Checking the quality of the reads
    4. MultiQC
2. Somatic short variant discovery analysis and annotation
    1. GATK re-processing (from fastq to rdy-to-use BAM : generate uBAM, MarkIlluminaAdapters, uBAM to fastq, BWA MEM, Merge BAM and UBAM, Mark Duplicates) - per sample.
    2. GATK_gvcf (HaplotypeCaller, from BAM to g.vcf ) - per sample per chromosome.
    3. GATK_merge_gvcf - CombineGVCFs combine g.vcf files to cohorts.
    4. GenotypeGVCFs - Perform joint genotyping on gVCF files produced by HaplotypeCaller, generate multi VCF file - per-cohort per-chromosome.
    5. GATK_hard_filters - Filter the multi VCF file - per cohort per chromosome.
    6. VEP - annotate the multi VCF file (`Variant Effect Predictor <https://www.ensembl.org/info/docs/tools/vep/index.html>`_. )
    7. GATK_SelectVariants_VEPfiltered - Separate multi VCF per-chromosome to one VCF per-samples per-chromosomes.
    8. GATK_CatVariants - Concatenate chromosome to get one VCF file for each sample.

        
Workflow Schema
~~~~~~~~~~~~~~~~

.. image:: GATK_workflow.jpg
   :alt: GATK workflow DAG

Requires
~~~~~~~~

`fastq` files. Paired end or single-end.

Programs required
~~~~~~~~~~~~~~~~~~

* `FastQC       <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
* `trimmomatic  <http://www.usadellab.org/cms/?page=trimmomatic>`_
* `GATK         <https://software.broadinstitute.org/gatk/>`_
* `VEP          <https://www.ensembl.org/info/docs/tools/vep/index.html>`_
* `Picard tools <https://broadinstitute.github.io/picard/>`_
* `BWA          <http://bio-bwa.sourceforge.net/>`_

Example of Sample File
~~~~~~~~~~~~~~~~~~~~~~

::

    Title	Metagenomics

    #SampleID	Type	Path    lane
    Sample1	Forward	/path/to/Sample1_F1.fastq.gz
    Sample1	Reverse	/path/to/Sample1_R1.fastq.gz
    Sample2	Forward	/path/to/Sample2_F1.fastq.gz
    Sample2	Reverse	/path/to/Sample2_R1.fastq.gz


Download
~~~~~~~~~

The workflow file is available :download:`here <../../../Workflows/GATK_workflow.yaml>`

