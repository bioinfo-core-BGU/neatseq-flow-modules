# -*- coding: UTF-8 -*-
""" 
``bwa_builder`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running bwa index builder:

Builds a bwa index for a fasta file stored at the project or sample level.

Determine which one will be used by specifying ``scope`` as either ``project`` or ``sample``.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output index files in one of the following slots:
    * ``self.sample_data[<sample>]["bwa_index"]``
    * ``self.sample_data["bwa_index"]``

* Puts the fasta file in one of the following slot:
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

    bwa_bld_ind:
        module: bwa_builder
        base: spades1
        script_path: /path/to/bwa index
        scope: project

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Li, H. and Durbin, R., 2009. **Fast and accurate short read alignment with Burrows–Wheeler transform**. *Bioinformatics*, 25(14), pp.1754-1760.
"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_bwa_builder(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        if not "scope" in self.params.keys():
            raise AssertionExcept("You must specify scope as 'project' or 'sample'\n")
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            # Initializing project bwa slot
            try:
                self.sample_data["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("Project does not have a nucl fasta defined. Check your 'scope'\n", sample)
            else:
                if "bwa_index" in self.sample_data.keys():
                    raise AssertionExcept("bwa index already seems to exist.\n")
            
                

        elif self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a nucl fasta defined. Can't build index\n", sample)
                else:
                    if "bwa_index" in self.sample_data[sample].keys():
                        raise AssertionExcept("bwa index already exists for sample.\n", sample)
            
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
            # self.sample_data["fasta.nucl"]
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
                output_prefix = sample + "_bwa_index"
                
                # Get constant part of script:
                self.script += self.get_script_const()

                # Add target for index
                self.script += "-p %s \\\n\t" % (use_dir + output_prefix)
                # Add source for index
                self.script += "%s \n\n" % self.sample_data[sample]["fasta.nucl"]


                self.sample_data[sample]["bwa_index"] = (sample_dir + output_prefix)
                self.sample_data[sample]["bwa_fasta"] = self.sample_data[sample]["fasta.nucl"]
                
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

                for ext in "amb ann bwt pac sa".split(" "):
                    self.stamp_file("%s%s.%s" % (sample_dir, output_prefix, ext))
                
                self.create_low_level_script()
                        
        
        else:  # scope == "project"
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            output_prefix = self.sample_data["Title"] + "_bwa_index"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            
            # Add target for index
            self.script += "-p %s \\\n\t" % (use_dir + output_prefix)
            # Add source for index
            self.script += "%s \n\n" % self.sample_data["fasta.nucl"]


            self.sample_data["bwa_index"] = (self.base_dir + output_prefix)
            self.sample_data["bwa_fasta"] = self.sample_data["fasta.nucl"]
            
            
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
       
        
            for ext in "amb ann bwt pac sa".split(" "):
                self.stamp_file("%s%s.%s" % (self.base_dir, output_prefix, ext))
            

            
            self.create_low_level_script()
                    

        
