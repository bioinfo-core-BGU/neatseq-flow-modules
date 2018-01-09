# -*- coding: UTF-8 -*-
""" 
``Trinity_gene_to_trans_map``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for running ``BUSCO``.




Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
Output:
~~~~~~~~~~~~~



Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_Trinity_gene_to_trans_map(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "BUSCO"
        
        
        if "scope" not in self.params:
            raise AssertionExcept("Please specify a 'scope': Either 'sample' or 'project'.")
        
                
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """


        if self.params["scope"] == "sample":
            # Check that "fasta" and "assembly" exist (signs that trinity has been executed)
            for sample in self.sample_data["samples"]:
                if ("fasta.nucl") not in self.sample_data[sample]:
                    raise AssertionExcept("It seems there is no sample-wide assembly file.", sample)
        elif self.params["scope"] == "project":
            # print self.sample_data.keys()
            if ("fasta.nucl") not in self.sample_data.keys():
                raise AssertionExcept("It seems there is no project-wide assembly file.")
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'.")


    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """

        pass

            
         

    def build_scripts(self):
    
        if self.params["scope"] == "project":
            self.build_scripts_project()
        else:
            self.build_scripts_sample()
            
            
    def build_scripts_project(self):
        
        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        self.script += self.get_script_const()

        output_bn = "%s.gene_trans_map" % os.path.basename(self.sample_data["fasta.nucl"])
        
        # The results will be put in data/step_name/name/Title
        self.script += "%s \\\n\t"  % self.sample_data["fasta.nucl"]
        self.script += "> %s \n\n" % os.path.join(use_dir,output_bn)
            
        

        # Store results to fasta and assembly slots:
        self.sample_data["gene_trans_map"] = os.path.join(self.base_dir,output_bn)
        self.stamp_file(self.sample_data["gene_trans_map"])
        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     

        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            

        
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            self.script += self.get_script_const()

            output_bn = "%s.gene_trans_map" % os.path.basename(self.sample_data[sample]["fasta.nucl"])
            
            # The results will be put in data/step_name/name/Title
            self.script += "%s \\\n\t"  % self.sample_data[sample]["fasta.nucl"]
            self.script += "> %s \n\n" % os.path.join(use_dir,output_bn)
            
        

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["gene_trans_map"] = os.path.join(self.base_dir,output_bn)
            self.stamp_file(self.sample_data[sample]["gene_trans_map"])
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
         

            self.create_low_level_script()
            
            
            
            
                 
            