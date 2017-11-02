# -*- coding: UTF-8 -*-
""" 
``IGV_toTDF`` (Included in main NeatSeq-Flow repo)
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running IGVtools toIGV:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* WIG file in the following slot:

    * ``sample_data[<sample>]["wig"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output tdf file in one the following slots:
    
    * ``self.sample_data[<sample>]["tdf"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "genome", "", "Path to chrom.sizes file for reference genome"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::

    IGV2TDF:
        module: IGV_toTDF
        base: samtools1
        script_path: /path/to/bin/java -Xmx1500m -jar /path/to/igvtools.jar toTDF 
        genome: /path/to/genome.chrom.sizes
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thorvaldsdóttir, H., Robinson, J.T. and Mesirov, J.P., 2013. **Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration**. *Briefings in bioinformatics*, 14(2), pp.178-192.        
"""

import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_IGV_toTDF(Step):

   
    
    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        # if "format" not in self.params.keys():
            # self.params["format"] = "tdf"
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Assert there is mapping data and a sorted bam in particular:
        for sample in self.sample_data["samples"]:      #Getting list of samples out of samples_hash
            if not "sam" in self.sample_data[sample] and not "bam" in self.sample_data[sample]:
                raise AssertionExcept("No BAM/SAM files exist\n", sample)
            if not "wig" in self.sample_data[sample]:
                raise AssertionExcept("No wig file exists", sample)
            # This might change if it is to be run on a snp file...

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
        
            
    
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            # self.spec_script_name = "_".join([self.step,self.name,sample])
            self.set_spec_script_name(sample)
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define input file (in future might work on snp files. this will have to be user-controlled?)
            input_file = self.sample_data[sample]["wig"]
            
            output_file = os.path.basename(input_file)

            # Get constant part of script:
            self.script += self.get_script_const()
            # Input file:
            self.script += "%s \\\n\t" % input_file
            # Output file:
            self.script += "%s%s.%s \\\n\t" % (use_dir, output_file, "tdf")

            # Add genome:
            self.script += "%s \n\n" % self.params["genome"]


            self.sample_data[sample]["tdf"] = "%s%s.%s" % (sample_dir, output_file, "tdf")
            self.stamp_file(self.sample_data[sample]["tdf"])
    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
        
