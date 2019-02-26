# -*- coding: UTF-8 -*-
"""
``GATK_pre_processing``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for generating *ready-to-GATK-use* BAM files from fastq files.

.. attention:: The module lacks the "base recalibration process (BQSR)" step

The programs included in the module are the following:

* ``FastqToSam`` Picard tool to generate uBAM
* ``MarkIlluminaAdapters`` Picard tool to Mark Illumina Adapters
* ``SamToFastq`` Picard tool uBAM to fastq
* ``MergeBamAlignment`` Picard tool to merge BAM and uBAM
* ``MarkDuplicates`` Picard tool to remove PCR duplicates
* ``BWA MEM`` mapping with BWA MEM



Requires
~~~~~~~~~~~~~~~

* A fastq file in the following locations:

    * ``self.sample_data[sample]["fastq.F"]``
    * ``self.sample_data[sample]["fastq.R"]``


Output
~~~~~~~~~~~~~~~~

    * ``self.sample_data[sample]["bam"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "picard_path", "path to PICARD", "Full path to the PICARD .jar file"
    "bwa_mem_path", "", ""
    "genome_reference", "", ""

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    GATK_pre_processing:
        module: GATK_pre_processing
        base: fQC_trim
        script_path: /path/to/java -jar /path/to/GenomeAnalysisTK.jar
        picard_path:     /path/to/picard.jar
        bwa_mem_path:    /path/to/bwa mem
        genome_reference:    /path/to/gatk/bundle/b37/human_g1k_v37_decoy.fasta
        threads: 20
        qsub_params:
            -pe: shared 20



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
http://broadinstitute.github.io/picard/

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

class Step_GATK_pre_processing(Step):
    """ A class that defines a pipeline step name (=instance).
    """


    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
            
            
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Prepare a list to store the qsub names of this steps scripts (will then be put in pipe_data and returned somehow)
        self.qsub_names=[]
        
       
        # Each iteration must define the following class variables:
            # spec_qsub_name
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample])
            self.script = ""
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            my_string = """
cd %(sample_dir)s
echo '\\n---------- generate uBAM -------------\\n'
%(picard_path)s FastqToSam \\
FASTQ=%(Forword_reads)s \\
FASTQ2=%(Revers_reads)s \\
OUTPUT=%(output_uBAM)s \\
READ_GROUP_NAME=%(sample_name)s \\
SAMPLE_NAME=%(sample_name)s \\
LIBRARY_NAME=%(group_name)s \\
PLATFORM_UNIT=Unit1 \\
PLATFORM=illumina \\
SEQUENCING_CENTER=BI \\
RUN_DATE=2016-08-20T00:00:00-0400


echo '\\n---------- MarkIlluminaAdapters -------------\\n'
%(picard_path)s MarkIlluminaAdapters \\
I=%(output_uBAM)s \\
O=%(output_markilluminaadapters)s \\
M= %(matrix_markilluminaadapters)s \\
TMP_DIR=%(sample_dir)s

echo '\\n---------- uBAM to fastq -------------\\n'
%(picard_path)s SamToFastq \
I=%(output_markilluminaadapters)s \\
FASTQ=%(output_samtofastq)s \\
CLIPPING_ATTRIBUTE=XT \\
CLIPPING_ACTION=2 \\
INTERLEAVE=true \\
NON_PF=true \\
TMP_DIR=%(sample_dir)s

echo '\\n---------- BWA MEM -------------\\n'
%(path_bwa_mem)s -M -t %(thread_number)s -p %(genome_reference)s \\
%(output_samtofastq)s > %(output_bwa_mem)s

echo '\\n---------- Merge BAM and UBAM -------------\\n'
%(picard_path)s MergeBamAlignment \\
R=%(genome_reference)s \\
UNMAPPED_BAM=%(output_uBAM)s \\
ALIGNED_BAM=%(output_bwa_mem)s \\
O=%(output_merge_bam_ubam)s \\
CREATE_INDEX=true \\
ADD_MATE_CIGAR=true \\
CLIP_ADAPTERS=false \\
CLIP_OVERLAPPING_READS=true \\
INCLUDE_SECONDARY_ALIGNMENTS=true \\
MAX_INSERTIONS_OR_DELETIONS=-1 \\
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
ATTRIBUTES_TO_RETAIN=XS \\
TMP_DIR=%(sample_dir)s


echo '\\n---------- Mark dup -------------\\n'
%(picard_path)s MarkDuplicates \\
INPUT=%(output_merge_bam_ubam)s \\
OUTPUT=%(output_duplicates)s \\
METRICS_FILE=%(matrix_duplicates)s \\
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
CREATE_INDEX=true \\
TMP_DIR=%(sample_dir)s


                    
            """ % { "sample_dir" : sample_dir,
                    "picard_path" : self.params["picard_path"],
                    "Forword_reads" : self.sample_data[sample]["fastq.F"],
                    "Revers_reads" : self.sample_data[sample]["fastq.R"],
                    "output_uBAM" : sample_dir + sample + "_fastqtosam.bam",
                    "output_markilluminaadapters": sample_dir + sample + "_markilluminaadapters.bam",
                    "matrix_markilluminaadapters": sample_dir + sample + "_markilluminaadapters_metrics.txt",
                    "group_name" : "OCD",
                    "sample_name" : sample,
                    "output_samtofastq" : sample_dir + sample + "_samtofastq_interleaved.fq",
                    "path_bwa_mem" : self.params["bwa_mem_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "output_bwa_mem" : sample_dir + sample + "_bwa_mem.bam",
                    "output_merge_bam_ubam" : sample_dir + sample + "_merge.bam",
                    "output_duplicates" : sample_dir + sample + "_duplicates.bam",
                    "matrix_duplicates" : sample_dir + sample + "_matrix_duplicates.txt",
                    "thread_number" : self.params["threads"]

            }
            
            my_string += """
rm -rf {merge_bam} \
    {merge_bai} \
    {bwa_mem_bam} \
    {samtofastq_interleaved} \
    {markilluminaadapters} \
    {fastqtosam}
    
            """.format(merge_bam=sample_dir + sample + "_merge.bam",
                       merge_bai=sample_dir + sample + "_merge.bai",
                       bwa_mem_bam=sample_dir + sample + "_bwa_mem.bam",
                       samtofastq_interleaved=sample_dir + sample + "_samtofastq_interleaved.fq",
                       markilluminaadapters=sample_dir + sample + "_markilluminaadapters.bam",
                       fastqtosam=sample_dir + sample + "_fastqtosam.bam"
                       )
            
            
            
            self.script += my_string
            #self.get_script_env_path()
            
            

            self.sample_data[sample]["bam"] = sample_dir + sample + "_duplicates.bam"
            self.stamp_file(self.sample_data[sample]["bam"])
            
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                    
