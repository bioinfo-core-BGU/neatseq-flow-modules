# -*- coding: UTF-8 -*-
""" 
``qiime_pick_otus``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's ``pick_otus.py``


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A fasta file in:

    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting OTU table in: 

    * ``self.sample_data["project_data"]["otu_table"]``
    

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_pick_otu_1:
        module: qiime_pick_otus
        base: q_chimera_usrch
        script_path: '{Vars.qiime_path}/pick_otus.py'
        setenv: {Vars.qiime_env}
    

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.



"""




import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_qiime_pick_otus(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_pick_otus"
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        
        try:
            self.sample_data["project_data"]["fasta.nucl"]
        except KeyError:
            raise AssertionExcept("fasta file does not exist.\n")

        
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
        for key in self.params["redir_params"].keys():
            self.script += "%s %s \\\n\t" % (key,self.params["redir_params"][key])
        self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
        # self.script += "-o %s \n\n" % self.base_dir
        self.script += "-o %s \n\n" % use_dir

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

        

        # Store location of demultiplexed folder
        # the OTU table is named like the input fasta, with the extension (.fna) replaced with _otus.txt
        self.sample_data["project_data"]["otu_table"] = self.base_dir + \
                    os.path.basename(os.path.splitext(self.sample_data["project_data"]["fasta.nucl"])[0]) + \
                    "_otus.txt";
        

        self.stamp_file(self.sample_data["project_data"]["otu_table"])
        
        self.create_low_level_script()
                    
            