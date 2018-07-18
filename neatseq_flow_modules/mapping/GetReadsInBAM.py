# -*- coding: UTF-8 -*-
""" 
``GetReadsInBAM``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running STAR mapper:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
* If ``scope`` is set (must come after ``STAR_builder`` module which populates the required slots):
    
    * STAR index directories in:

        * ``sample_data[<sample>]["STAR_index"]``  if ``scope`` = "sample"
        * ``sample_data["STAR_index"]``  if ``scope`` = "project"

    * Reference fasta files in:

        * ``sample_data[<sample>]["STAR_fasta"]``  if ``scope`` = "sample"
        * ``sample_data["STAR_fasta"]``  if ``scope`` = "project"

    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts output sam files in the following slots:

    * ``self.sample_data[<sample>]["sam"]``

* Alternatively, if ``--outSAMtype`` is set to ``BAM``, puts output BAM files in the following slots:

    * ``self.sample_data[<sample>]["bam"]``
    * ``self.sample_data[<sample>]["bam_unsorted"]``

* High confidence collapsed splice junctions (SJ.out.tab  file) will be stored in:

    * ``self.sample_data[<sample>]["SJ.out.tab"]``

* If ``--quantMode`` contains ``TranscriptomeSAM``, alignments BAM translated into transcript coordinates will be stored in:

    * ``self.sample_data[<sample>]["TranscriptomeSAM"]``

* If ``--quantMode`` contains ``GeneCounts``, the ``ReadsPerGene.out.tab`` file will be stored:

    * ``self.sample_data[<sample>]["GeneCounts"]``

* If ``--outWigType`` is set, will store outputs in:

    * if ``--outWigType`` is ``wiggle``

        * ``self.sample_data[<sample>]["wig2_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["wig2_Unique"]``
        * ``self.sample_data[<sample>]["wig1_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["wig1_Unique"]``
        * ``self.sample_data[<sample>]["wig"]``

    * if ``--outWigType`` is ``bedGraph``

        * ``self.sample_data[<sample>]["bdg2_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["bdg2_Unique"]``
        * ``self.sample_data[<sample>]["bdg1_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["bdg1_Unique"]``
        * ``self.sample_data[<sample>]["bdg"]``
    
                
* Puts the name of the mapper in:
    ``self.sample_data[<sample>]["mapper"]``

* Puts fasta of reference genome (if one is given in param file) in:
    ``self.sample_data[<sample>]["reference"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ref_genome", "path to genome fasta", ""
    "scope", "project | sample", "The scope from which to take the genome directory"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

::

    STAR_map:
        module:             STAR_mapper
        base:               STAR_bld_ind
        script_path:        /path/to/STAR
        redirects:
            --readMapNumber:    1000
            --genomeDir:        /path/to/genome/STAR_index/
            
            
**For project STAR index:**

::

    STAR_map:
        module:             STAR_mapper
        base:               STAR_bld_ind
        script_path:        /path/to/STAR
        scope:              project
        redirects:
            --readMapNumber:    1000
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*, 29(1), pp.15-21.

"""


import os, re
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Menachem Sklarz"
__version__ = "1.2.0"

class Step_GetReadsInBAM(Step):

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        # self.auto_redirs = ["--readFilesCommand", "--readFilesIn", "--outFileNamePrefix", "--outTmpDir", "--outStd"]


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
        
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
 
            # Define location and prefix for output files:
            output_prefix = sample + "_mapped_reads.txt"
            self.script = ""
            

            self.script =  """
{script_path} view {bam} | cut -f 1 | sort | uniq > {reads}
""".format(script_path=self.params["script_path"],
           bam=self.sample_data[sample]["bam"], 
           reads=use_dir + output_prefix)

            # Storing name of mapper. might be useful:
            self.sample_data[sample]["read_list"] = sample_dir + output_prefix
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
