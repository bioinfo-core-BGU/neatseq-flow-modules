# -*- coding: UTF-8 -*-
""" 
``merge_table``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging sample tables into a single project wide table. The table can be with or without a header line. 

Can be used for merging fasta and fastq files as well.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A table file in any slot:

    * ``sample_data[<sample>][<file.type>]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slot:
        
    * ``sample_data[<file.type>]``




Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "", "A file type that exists in all samples."
    "script_path", "", "Leave blank"
    "header","","The number of header lines each table has. The header will be used for the complete table and all other headers will be removed. If there is no header line, set to 0 or leave out completely. **If not specified but not set, will default to 1!**"



Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    merge_blast_tables:
        module:         merge_table
        base:           merge1
        script_path:    
        type:           [blast,blast.prot]
        header:         1


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_merge_table(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        if "header" not in self.params:
            self.params["header"] = 0
        elif not self.params["header"]:
            self.params["header"] = 1
        else:
            pass
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        type = set()
        
        if "type" not in self.params:
            raise AssertionExcept("You must supply a file type or list thereof to merge")
        else: # If passed by user, converting even single values to a list.
            if not isinstance(self.params["type"], list):
                self.params["type"] = [self.params["type"]]

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            for type in self.params["type"]:
                if type not in self.sample_data[sample]:
                    raise AssertionExcept("Type %s does not exist in sample" % type, sample)

        
        
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
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"],type])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Define location and prefix for output files:
            output_fn = "{filename}".format(filename = ".".join([self.sample_data["Title"],type]))
            
            
            if self.params["header"]==0:
                self.script += "cat \\\n\t".format(header = self.params["header"])
            else:
                self.script += "sed -s '1,{header}d' \\\n\t".format(header = self.params["header"])

            
            # # Get constant part of script:
            # self.script += self.get_script_const()
            # # Files to merge:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s \\\n\t" % self.sample_data[sample][type]
            
            self.script += "> {dir}{file}\n\n".format(dir=use_dir,file=output_fn) 
            
          
            self.sample_data[type] = "%s%s" % (self.base_dir, output_fn)
            self.stamp_file(self.sample_data[type])
                    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                            
