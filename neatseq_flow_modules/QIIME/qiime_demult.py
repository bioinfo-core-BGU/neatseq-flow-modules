# -*- coding: UTF-8 -*-
"""
``qiime_demult``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's multiple_split_libraries_fastq.py:

The reads from step qiime_prep are combined into one seqs.fna file.

.. Note:: The module has not been tested on other types of data, such as undemultiplexed reads. It should work but there will probably be unexpected problems.



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A directory of read files with smaple names coded in the file names, such as the directory produced by qiime_prep:

    * ``sample_data["qiime.prep_links_dir"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting ``seqs.fna`` file in the following slots:

    * ``self.sample_data["project_data"]["qiime.demult_seqs"]``
    * ``self.sample_data["project_data"]["qiime.fasta"]``
    * ``self.sample_data["project_data"]["fasta.nucl"]``

    
    

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_demult_1:
        module: qiime_demult
        base: q_prep_1
        script_path: '/path/to/multiple_split_libraries_fastq.py'
        redirects:
            --demultiplexing_method: sampleid_by_file
            --include_input_dir_path: null
            --parameter_fp: /path/to/qiime_params
            --remove_filepath_in_name: null
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.



"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_qiime_demult(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_demult"
        
                
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # If does not exist 
        try:
            self.sample_data["project_data"]["qiime.prep_links_dir"]
        except KeyError:
            raise AssertionExcept("It seems like qiime_demult is the first qiime step. At the moment, it must come after qiime_prep...\n" )
        
        # # Testing the existance of a legitimate parameter file in the pipeline parameter file:
        # if not "--parameter_fp" in self.params["redir_params"].keys() and not "-p" in self.params["redir_params"].keys():
            # if not "qiime.parameters" in self.sample_data:
                # raise AssertionExcept("You must supply a parameter file\n")
            # else:
                # self.write_warning("Using parameter file passed previously")
        # else:
            # if "qiime.parameters" in self.sample_data:
                # self.write_warning("Overriding parameter file passed previously")
            # self.params["qiime.parameters"] = self.params["redir_params"]["--parameter_fp"] if "--parameter_fp" in self.params["redir_params"].keys() else self.params["redir_params"]["-p"]

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      

    def build_scripts(self):
        

        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        
        # if "local" in self.params.keys():
            # local_dir = "/local/bioinfo/" + "_".join([self.step,self.name,self.pipe_data["run_code"]]) + os.sep
            # self.script += "mkdir -p %s \n\n" % local_dir
            # use_dir = local_dir 
        # else:
            # use_dir = self.base_dir

        ### Step 1b: Adding demultiplexing tyo script:
        self.script += self.get_script_const()        # Gets the "env", "script_path" and "redir_params" part of the script which is always the same...
        self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["qiime.prep_links_dir"]
        # self.script += "-o %s \n\n" % self.base_dir
        self.script += "-o %s \n\n" % use_dir


        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


        # Store location of demultiplexed folder
        #   Static: Stores the position of the demulted fasta throughout the pipeline
        self.sample_data["project_data"]["qiime.demult_seqs"] = self.base_dir + "seqs.fna";
        #   Dynamic: Later steps that modify the fasta will set this to the modified filename
        # self.sample_data["project_data"]["qiime.fasta"] = self.base_dir + "seqs.fna";
        self.sample_data["project_data"]["fasta.nucl"]  = self.base_dir + "seqs.fna";
        

        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            