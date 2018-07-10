# -*- coding: UTF-8 -*-
""" 
``TransDecoder`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``TransDecoder`` on a transcripts file.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
``fasta`` files in at least one of the following slots:
    
    * ``sample_data[<sample>]["fasta.nucl"]``  (if ``scope = sample``)
    * ``sample_data["fasta.nucl"]``  (if ``scope = project``)

    
Output:
~~~~~~~~~~~~~

* If ``scope = project``:

    * Protein fasta in ``self.sample_data["fasta.prot"]``
    * Gene fasta in ``self.sample_data["fasta.nucl"]``
    * Original transcripts in ``self.sample_data["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data["gff3"]``

* If ``scope = sample``:

    * Protein fasta in ``self.sample_data[<sample>]["fasta.prot"]``
    * Gene fasta in ``self.sample_data[<sample>]["fasta.nucl"]``
    * Original transcripts in ``self.sample_data[<sample>]["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data[<sample>]["gff3"]``


                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Determine weather to use sample or project transcripts file."
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_Transdecode_highExpr:
        module:             TransDecoder
        base:               Split_Fasta
        script_path:        {Vars.paths.TransDecoder}
        scope:              sample
        
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_TransDecoder(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".Trinity.fasta"
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                if "fasta.nucl" not in self.sample_data:
                    raise AssertionExcept("Project does not have a nucl fasta.")

            elif self.params["scope"]=="sample":
                
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    if "fasta.nucl" not in self.sample_data[sample]:
                        raise AssertionExcept("Sample does not have a nucl fasta.", sample)
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
        
         

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
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

        output_dir_base = "%s.transdecoder_dir" % os.path.basename(self.sample_data["fasta.nucl"])
         
        self.script = "cd {usedir}\n\n".format(usedir=use_dir)
        self.script += self.get_script_const()
        self.script += "-t {fasta_nucl}\n\n".format(fasta_nucl = self.sample_data["fasta.nucl"])
        
        self.script += "cd {home_dir}\n\n".format(home_dir = self.pipe_data["home_dir"])

        

        # Store results to fasta and assembly slots:
        self.sample_data["fasta.prot"] = os.path.join(self.base_dir, output_dir_base, "longest_orfs.pep")
        self.sample_data["transcripts.fasta.nucl"] = self.sample_data["fasta.nucl"]
        self.sample_data["fasta.nucl"] = os.path.join(self.base_dir, output_dir_base, "longest_orfs.cds")
        self.sample_data["gff3"] = os.path.join(self.base_dir, output_dir_base, "longest_orfs.gff3")
        
        self.stamp_file(self.sample_data["fasta.prot"])
        self.stamp_file(self.sample_data["fasta.nucl"])
        self.stamp_file(self.sample_data["gff3"])

        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     
            
        
        
        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            output_dir_base = "%s.transdecoder_dir" % os.path.basename(self.sample_data[sample]["fasta.nucl"])

            self.script = "cd {usedir}\n\n".format(usedir=use_dir)
            self.script += self.get_script_const()
            self.script += "-t {fasta_nucl}\n\n".format(fasta_nucl = self.sample_data[sample]["fasta.nucl"])
            
            self.script += "cd {home_dir}\n\n".format(home_dir = self.pipe_data["home_dir"])

            

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["fasta.prot"] = os.path.join(sample_dir, output_dir_base, "longest_orfs.pep")
            self.sample_data[sample]["transcripts.fasta.nucl"] = self.sample_data[sample]["fasta.nucl"]
            self.sample_data[sample]["fasta.nucl"] = os.path.join(sample_dir, output_dir_base, "longest_orfs.cds")
            self.sample_data[sample]["gff3"] = os.path.join(sample_dir, output_dir_base, "longest_orfs.gff3")
            
            self.stamp_file(self.sample_data[sample]["fasta.prot"])
            self.stamp_file(self.sample_data[sample]["fasta.nucl"])
            self.stamp_file(self.sample_data[sample]["gff3"])

            
                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                        
            
            
                 
            
     