# -*- coding: UTF-8 -*-
""" 
``genomeCoverageBed`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for running bedtools genomecov:

The module builds a bedgraph (bdg) file based on an existing BAM file.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BAM file in the following slot:

    * ``sample_data[<sample>]["bam"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output BedGraph files in the following slots:
    * ``sample_data[<sample>]["bdg"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "-g", "path to chrom.sizes", "You must redirect the -g parameter. Create the chrom.sizes file for the reference genome with ``samtools faidx`` followed by ``cut -f1,2``."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    genCovBed_bwt1:
        module: genomeCoverageBed
        base: sam_bwt1
        script_path: /path/to/bedtools/bin/genomeCoverageBed
        redirects:
            -bg: 
            -g: /path/to/ref_genome/ref_genome.chrom.sizes



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
    
    

import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_genomeCoverageBed(Step):
  
    
    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        if not "-g" in self.params["redir_params"]:
            raise AssertionExcept("You must pass a genome (-g) parameter with chromosome sizes!!!")
            
        if "-bg" not in self.params["redir_params"]:
            self.params["redir_params"]["-bg"] = None
            print "In step %s: Adding -bg to genomeCoverageBed parameters!!!" % self.name
 
        if "-ibam" in self.params["redir_params"]:
            raise AssertionExcept("Please DO NOT pass -ibam. That is done automatically")
            

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Assert there is mapping data and a sorted bam in particular:
        for sample in self.sample_data["samples"]:      #Getting list of samples out of samples_hash
            
            if not "bam" in self.sample_data[sample].keys():
                raise AssertionExcept("No BAM file defined.\nDid you run samtools first?\n", sample)

            
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
        
        
        # Each iteration must define the following class variables:
            # self.spec_script_name
            # self.script
        
        # self.base_dir    
    
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define input file
            input_file = self.sample_data[sample]["bam"]
            
            output_file = "%s.bdg" % os.path.basename(input_file)

            # SCript path to genomeCoverageBed
            self.script += self.get_script_const()
            self.script += "-ibam %s > \\\n\t" % input_file
            self.script += "%s\n\n" % (use_dir + output_file)
            
            
            self.sample_data[sample]["bdg"] = "%s%s" % (sample_dir, output_file)
            self.stamp_file(self.sample_data[sample]["bdg"])
            
    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
        
