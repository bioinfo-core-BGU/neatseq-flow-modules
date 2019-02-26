# -*- coding: UTF-8 -*-
""" 
``qiime_align_seqs``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``QIIME's align_seqs.py``:

Can be used for the parallel versions thereof: ``parallel_align_seqs_pynast.py``

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A fasta file in:

    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting aligned fasta file in:

    * ``self.sample_data["project_data"]["fasta.nucl"]``
    * ``self.sample_data["project_data"]["fasta.aligned"]``
    
* Stores the old, unaligned version in:

    * ``self.sample_data["project_data"]["fasta.unaligned"]``


    
.. Note:: When using ``parallel_align_seqs_pynast.py``, the module tries to build the scripts appropriately. **It is wise to check the parallel scripts before running them...**
  
  


    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_align_para:
        module: qiime_align_seqs
        base: q_rep_set_1
        script_path: '{Vars.qiime_path}/parallel_align_seqs_pynast.py'
        setenv: {Vars.qiime_env}
        redirects:
            --jobs_to_start: 5
            --retain_temp_files: 
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.



"""




import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_qiime_align_seqs(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_align_seqs"
            
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["project_data"]["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like qiime_demult is the first qiime step. At the moment, it must come after qiime_prep...\n" )

        try:
            self.sample_data["project_data"]["fasta.nucl"]
        except KeyError:
            raise AssertionExcept("fasta file does not exist.\n")

        try:
            self.sample_data["project_data"]["otu_table"]
        except KeyError:
            self.write_warning("otu table does not exist.\n")


        
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


        ### Step 1b: Adding demultiplexing tyo script:
        # if "env" in self.params.keys():         # Add optional environmental variables.
            # self.script += "env %s \\\n\t" % self.params["env"]
        # self.script += self.params["script_path"] + " \\\n\t"
        self.script += self.get_script_env_path()
        if re.search("parallel",self.params["script_path"]):
            # Adding --poll_directly: This will keep the original job "alive" while subprocesses are running. Important for keeping the downstream steps waiting for completion
            self.script += "--job_prefix %s \\\n\t--poll_directly \\\n\t"  % self.spec_script_name

        for key in list(self.params["redir_params"].keys()):
            self.script += "%s %s \\\n\t" % (key,self.params["redir_params"][key] if self.params["redir_params"][key] else "")

        if not set(["-t","--template_fp"]) & set(self.params["redir_params"].keys()):      # no template_fp was passed in redir_params:
            self.write_warning("Template file not defined. Using default (find it with print_qiime_config.py\n")
    
        self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
        self.script += "-o %s \n\n" % use_dir


        outfile = os.path.basename(self.sample_data["project_data"]["fasta.nucl"])
        outfile = re.sub("\.(fas|fasta|fna|fa)$","",outfile)

        # Store location of demultiplexed folder
        self.sample_data["project_data"]["fasta.aligned"] = self.base_dir + outfile + "_aligned.fasta";    # Saving aligned fasta separately
        self.sample_data["project_data"]["fasta.unaligned"] = self.sample_data["project_data"]["fasta.nucl"]  # Saving unaligned fasta in case someone wants it...
        self.sample_data["project_data"]["fasta.nucl"] = self.sample_data["project_data"]["fasta.aligned"]      # Using aligned as active fasta file.

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

        self.stamp_file(self.sample_data["project_data"]["fasta.aligned"])
        
        self.create_low_level_script()
                    
            