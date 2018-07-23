# -*- coding: UTF-8 -*-
""" 
``GetReadsInBAM``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for extracting read names from a bam file.

.. Tip:: Can be used together with ``FilterSamReads`` from the PICARD tools to extract a subset of the reads.
   The output in slot ``read_list`` can be used in the ``FilterSamReads`` ``READ_LIST_FILE`` parameter.
   At the moment, you will have to do this with a Generic instance.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Sample scope SAM or BAM files in one of the following slots (if both exist, BAM will be used):

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["sam"]``
    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output reads list in the following slots:

    * ``self.sample_data[<sample>]["read_list"]``


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

        for sample in self.sample_data["samples"]:
            if "bam" not in self.sample_data[sample] and "sam" not in self.sample_data[sample]:
                raise AssertionExcept("No BAM or SAM files defined for sample", sample)


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

            if "bam" in self.sample_data[sample]:
                input_file = self.sample_data[sample]["bam"]
            else:
                input_file = sample_data[sample]["sam"]

            self.script =  """
{script_path} view {bam} | cut -f 1 | sort | uniq > {reads}
""".format(script_path=self.params["script_path"],
           bam=input_file,
           reads=use_dir + output_prefix)

            # Storing name of mapper. might be useful:
            self.sample_data[sample]["read_list"] = sample_dir + output_prefix
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
