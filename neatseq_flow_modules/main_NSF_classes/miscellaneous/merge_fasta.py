# -*- coding: UTF-8 -*-
""" 
``merge_fasta``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Attention:: This module replaces the ``merge_project`` module. It will be one of a collection of ``merge_XXX`` modules.


A module for merging sample `fasta` files into a single project wide `fasta` file.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A `fasta` file in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slots:
        
    * ``sample_data["fasta.nucl"]``
    * ``sample_data["fasta.prot"]``



Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "nucl|prot", "Type of fasta to merge. If ommitted, both will be merged in parallel (if they exist)."
    "script_path", "cat", "This must be *cat*. Do not change it..."



Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    merge_proj:
        module:         merge_fasta
        base:           merge1
        script_path:    cat
        type:           fasta.nucl


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"

class Step_merge_fasta(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        type = set()
        
        if "type" not in self.params:
            # self.params["type"] = ["fastq.S","fastq.F","fastq.R","nucl","prot"]
            self.params["type"] = ["fasta.nucl","fasta.prot"]
        else: # If passed by user, converting even single values to a list.
            if not isinstance(self.params["type"], list):
                self.params["type"] = [self.params["type"]]

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            self.params["type"] = [x for x in self.params["type"] if x in list(self.sample_data[sample].keys())]
            if not self.params["type"]:
                raise AssertionExcept("None of the input types exist in sample", sample)

        
        
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

        for type in self.params["type"]:
            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],type])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Define location and prefix for output files:
            ext_dict = {"fasta.nucl":"fna", "fasta.prot":"faa", "fastq.S":"fq", "fastq.F":"fq", "fastq.R":"fq"}
            output_fn = "{filename}.{ext}".format(filename = "_".join([self.sample_data["Title"],type]), \
                                                 ext = ext_dict[type])
            
            
            # Get constant part of script:
            self.script += self.get_script_const()
            # Files to merge:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s \\\n\t" % self.sample_data[sample][type]
            
            self.script += "> %s%s" % (use_dir, output_fn)
            
          
            self.sample_data["project_data"][type] = "%s%s" % (self.base_dir, output_fn)
            self.stamp_file(self.sample_data["project_data"][type])
                    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                            
