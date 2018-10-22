# -*- coding: UTF-8 -*-
""" 
``bowtie2_builder`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running bowtie2 index builder:

Builds a bowtie2 index for a fasta file stored at the project or sample level.

Determine which one will be used by specifying ``scope`` as either ``project`` or ``sample``.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output index files in one of the following slots:
    * ``self.sample_data[<sample>]["bowtie2_index"]``
    * ``self.sample_data["project_data"]["bowtie2_index"]``

* Puts the fasta file in the following slot:
    * ``self.sample_data[<sample>]["reference"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Indicates whether to use a project fasta or a sample fasta."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    bwt2_build:
        module: bowtie2_builder
        base: trinity1
        script_path: /path/to/bowtie2-build
        scope: project

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Langmead, B. and Salzberg, S.L., 2012. **Fast gapped-read alignment with Bowtie 2**. *Nature methods*, 9(4), pp.357-359.

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_bowtie2_builder(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        if not "scope" in self.params.keys():
            raise AssertionExcept("You must specify scope as 'project' or 'sample'\n")
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            # Initializing project bowtie2 slot
            try:
                self.sample_data["project_data"]["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("Project does not have a nucl fasta defined. Check your 'scope'\n", sample)
            # else:
                # if "bowtie2_index" in self.sample_data.keys():
                    # raise AssertionExcept("bowtie2 index already seems to exist.\n")
            
                

        elif self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a nucl fasta defined. Can't build index\n", sample)
                else:
                    if "bowtie2_index" in self.sample_data[sample].keys():
                        raise AssertionExcept("bowtie2 index already exists for sample.\n", sample)
            
        else:
            raise AssertionExcept("Scope must be either 'sample' or 'project'")
                
             
        
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
                output_prefix = use_dir + sample + "_bowtie2_index"
                
                # Get constant part of script:
                self.script += self.get_script_const()
                
                self.script += "%s \\\n\t" % self.sample_data[sample]["fasta.nucl"]
                self.script += "%s \n\n" % output_prefix


                self.sample_data[sample]["bowtie2_index"] = output_prefix
                self.sample_data[sample]["bowtie2_fasta"] = self.sample_data[sample]["fasta.nucl"]
                # self.stamp_dir_files(sample_dir)
        
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
           
                
                
                self.create_low_level_script()
                        
        
        else:  # scope == "project"
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()

            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            output_prefix = use_dir + self.sample_data["Title"] + "_bowtie2_index"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            
            self.script += "%s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
            self.script += "%s \n\n" % output_prefix


            self.sample_data["project_data"]["bowtie2_index"] = output_prefix
            self.sample_data["project_data"]["bowtie2_fasta"] = self.sample_data["project_data"]["fasta.nucl"]

        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    

        
