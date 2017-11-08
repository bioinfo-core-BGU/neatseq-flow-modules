Preparation and QC
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: main_NSF_classes.preparing.merge

.. automodule:: main_NSF_classes.preparing.fastqc_html

.. automodule:: main_NSF_classes.preparing.trimmo

.. automodule:: main_NSF_classes.Reports.Multiqc

.. automodule:: Liron.Cutadapt_module.Cutadapt

.. automodule:: Liron.Trim_Galore_module.Trim_Galore


Mapping
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: main_NSF_classes.mapping.bowtie2_builder

.. automodule:: main_NSF_classes.mapping.bowtie2_mapper

.. automodule:: main_NSF_classes.mapping.bowtie1_builder

.. automodule:: main_NSF_classes.mapping.bowtie1_mapper

.. automodule:: main_NSF_classes.mapping.bwa_builder

.. automodule:: main_NSF_classes.mapping.bwa_mapper

.. automodule:: mapping.STAR_mapper

.. automodule:: mapping.STAR_builder

.. automodule:: main_NSF_classes.mapping.samtools

.. automodule:: main_NSF_classes.Reports.Multiqc




Genome Assembly
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1

.. automodule:: Assembly.clc_assembl

.. automodule:: Assembly.megahit_assembl

.. automodule:: main_NSF_classes.Assembly.spades_assembl


RNA-seq Aassembly
========================================

.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: main_NSF_classes.RNA_seq.trinity

.. automodule:: main_NSF_classes.RNA_seq.add_trinity_tags

.. automodule:: Liron.RSEM_module.RSEM

 
Assembly QC
========================================

.. contents:: Modules included in this section
   :local:
   :depth: 1

.. automodule:: main_NSF_classes.Assembly.quast

Modules for sequence annotation
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1

.. automodule:: Annotation.prokka_old

.. automodule:: Liron.Prokka_module.Prokka


Modules for sequence clustering
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1



.. automodule:: clustering.cd_hit

.. automodule:: clustering.vsearch_cluster

.. automodule:: clustering.vsearch_derepel

Modules for metagenomics
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1



.. automodule:: metagenomics.HUMAnN2

.. automodule:: metagenomics.kaiju

.. automodule:: metagenomics.kraken

.. automodule:: metagenomics.metaphlan2


Modules for microbiology
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1



.. automodule:: microbiology.CARD_RGI

.. automodule:: Liron.cgMLST_and_MLST_typing_module.cgMLST_and_MLST_typing

.. automodule:: Liron.Roary_module.Roary

.. automodule:: Liron.Snippy_module.Snippy





Modules for QIIME (version 1.9)
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: QIIME.qiime_prep

.. automodule:: QIIME.qiime_demult

.. automodule:: QIIME.qiime_chimera

.. automodule:: QIIME.qiime_pick_otus

.. automodule:: QIIME.qiime_pick_rep_set

.. automodule:: QIIME.qiime_align_seqs

.. automodule:: QIIME.qiime_filter_alignment

.. automodule:: QIIME.qiime_assign_taxonomy

.. automodule:: QIIME.qiime_make_phylogeny

.. automodule:: QIIME.qiime_make_otu_table

.. automodule:: QIIME.qiime_filter_samples_from_otu_table

.. automodule:: QIIME.qiime_filter_otus

.. automodule:: QIIME.qiime_sort_otu_table

.. automodule:: QIIME.qiime_divers



Modules for various reporting programs
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1



.. automodule:: Reports.NGSplot

.. automodule:: main_NSF_classes.Reports.Multiqc

.. automodule:: Liron.Collect_results_module.Collect_results




Modules for sequence-searching related tasks
============================================



.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: main_NSF_classes.searching.makeblastdb

.. automodule:: main_NSF_classes.searching.blast

.. automodule:: searching.blast_new

.. automodule:: searching.parse_blast

.. automodule:: Liron.Gassst_module.Gassst

BAM conversion to other formats
============================================

.. contents:: Modules included in this section
   :local:
   :depth: 1

.. automodule:: main_NSF_classes.mapping.genomeCoverageBed

.. automodule:: main_NSF_classes.UCSCtools.UCSC_BW_wig

.. automodule:: main_NSF_classes.IGVtools.IGV_count

.. automodule:: main_NSF_classes.IGVtools.IGV_toTDF



Modules for various variant-related tasks
============================================



.. contents:: Modules included in this section
   :local:
   :depth: 1



.. automodule:: variants.freebayes

.. automodule:: variants.mpileup_varscan

.. automodule:: main_NSF_classes.variants.vcftools


Modules for ChIP-seq
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1

.. automodule:: main_NSF_classes.ChIP_seq.macs2_callpeak

.. automodule:: ChIP_seq.macs2_bdgcmp

.. automodule:: ChIP_seq.CEAS

Generic module
========================================


.. contents:: Modules included in this section
   :local:
   :depth: 1


.. automodule:: main_NSF_classes.Generic.Generic

