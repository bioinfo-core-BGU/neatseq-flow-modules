# -*- coding: UTF-8 -*-
""" 
``STAR_builder``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running STAR genome index construction:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.nucl"]``
       

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Puts output index files in one of the following slot:

    * ``self.sample_data[<sample>]["STAR_index"]``
    * ``self.sample_data["STAR_index"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Not used"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::


    STAR_bld_ind:
        module:             STAR_builder
        base:               trinity1
        script_path:        /path/to/STAR
        scope:              project
        qsub_params:
            queue:          star.q
        redirects:
            --genomeSAindexNbases:  12
            --genomeChrBinNbits:    10

    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*, 29(1), pp.15-21.

"""


import os, re
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_STAR_builder(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"


        
        for redir2remove in ["--runMode", "--genomeDir", "--genomeFastaFiles"]:
            if redir2remove in self.params["redir_params"]:
                del self.params["redir_params"][redir2remove]
                self.write_warning("You are not supposed to specify %s in redirects. We set it automatically" % redir2remove)


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" not in self.params.keys():
            # Try guessing scope:
            try:  # Does a nucl fasta exist for project?
                self.sample_data["fasta.nucl"]
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
        if self.params["scope"] == "sample":
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
     
                # Get constant part of script:
                self.script += self.get_script_const()
        
                self.script += "--runMode genomeGenerate \\\n\t"
                self.script += "--genomeDir %s \\\n\t"  % use_dir
                self.script += "--genomeFastaFiles %s \n\n"  % self.sample_data[sample]["fasta.nucl"]


                self.sample_data[sample]["STAR_index"] = sample_dir
                self.sample_data[sample]["STAR_fasta"] = self.sample_data[sample]["fasta.nucl"]
        
            
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
           
                
                
                self.create_low_level_script()
        else:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Get constant part of script:
            self.script += self.get_script_const()
            self.script += "--runMode genomeGenerate \\\n\t"
            self.script += "--genomeDir %s \\\n\t"  % use_dir
            self.script += "--genomeFastaFiles %s \n\n"  % self.sample_data["fasta.nucl"]



            self.sample_data["STAR_index"] = self.base_dir
            self.sample_data["STAR_fasta"] = self.sample_data["fasta.nucl"]
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    

        
        
