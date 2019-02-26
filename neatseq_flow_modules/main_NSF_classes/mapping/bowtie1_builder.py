# -*- coding: UTF-8 -*-
""" 
``bowtie1_builder`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running bowtie1 index builder:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.nucl"]``
    

output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
Puts output index files in one of the following slot:
    * ``self.sample_data[<sample>]["bowtie1.index"]``
    * ``self.sample_data["project_data"]["bowtie1.index"]``
            

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "path to bowtie1 index", "If not given, will look for a project bowtie1 index and then for a sample bowtie1 index"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    bwt1_bld_ind:
        module: bowtie1_builder
        base: trinity1
        script_path: /path/to/bowtie
        scope: project

    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Langmead, B., Trapnell, C., Pop, M. and Salzberg, S.L., 2009. **Ultrafast and memory-efficient alignment of short DNA sequences to the human genome**. *Genome biology*, 10(3), p.R25.

"""
    
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_bowtie1_builder(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        
        if "scope" not in list(self.params.keys()):
            # Try guessing scope:
            try:  # Does a nucl fasta exist for project?
                self.sample_data["project_data"]["fasta.nucl"]
            except KeyError:
                self.params["scope"] = "sample"
            else:
                self.params["scope"] = "project"
        else:
            # Check scope is legitimate
            if not self.params["scope"] in ["project","sample"]:
                raise AssertionExcept("Scope must be either 'sample' or 'project'\n")

        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a nucl fasta defined. Can't build index\n", sample)





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
        
        # try:    # Check if fasta nucl exists:
            # self.sample_data["project_data"]["fasta.nucl"]
            
        # except KeyError:   # If not, search in samples
        
        if self.params["scope"] == "sample":
        
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
 
                # Define location and prefix for output files:
                output_prefix = use_dir + sample + "_bowtie_index"
                
                # Get constant part of script:
                self.script += self.get_script_const()
                
                self.script += "%s \\\n\t" % self.sample_data[sample]["fasta.nucl"]
                self.script += "%s \n\n" % output_prefix


                self.sample_data[sample]["bowtie1.index"] = output_prefix
                self.sample_data[sample]["bowtie1.fasta"] = self.sample_data[sample]["fasta.nucl"]
                # self.stamp_dir_files(sample_dir)
        
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
           
                
                
                self.create_low_level_script()
                        
        
        else:  # If found, build on project fasta nucl
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()

            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            output_prefix = use_dir + self.sample_data["Title"] + "_bowtie_index"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            
            self.script += "%s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
            self.script += "%s \n\n" % output_prefix


            self.sample_data["project_data"]["bowtie1.index"] = output_prefix
            self.sample_data["project_data"]["bowtie1.fasta"] = self.sample_data["project_data"]["fasta.nucl"]

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    

        
