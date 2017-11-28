# -*- coding: UTF-8 -*-
""" 

``qiime_make_phylogeny``
--------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's ``make_phylogeny.py``


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A fasta file in:

    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting OTU table in: 

    * ``self.sample_data["phylotree"]``
    
    

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_phylo_1:
        module: qiime_make_phylogeny
        base: q_filt_align_1
        script_path: '{Vars.qiime_path}/make_phylogeny.py'
        setenv: {Vars.qiime_env}
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.



"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_qiime_make_phylogeny(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_make_phylogeny"
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like this is the first qiime step. At the moment, it must come after qiime_prep...\n" )
        
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      

    def build_scripts(self):
        

        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        if "fasta.aligned" not in self.sample_data.keys():
            raise AssertionExcept("You are trying to run 'make_phylogeny' on an unaligned fasta file!\n")
        outfile = os.path.basename(self.sample_data["fasta.nucl"])
        outfile = re.sub("\.(fas|fasta|fna|fa)$","",outfile) + ".tre"
        logfile = ".".join([outfile,"log"])

        ### Step 1b: Adding demultiplexing tyo script:
        self.script += self.get_script_const()        # Gets the "env", "script_path" and "redir_params" part of the script which is always the same...

        self.script += "-i %s \\\n\t" % self.sample_data["fasta.nucl"]
        # self.script += "-o %s \n\n" % self.base_dir
        self.script += "-o %s \\\n\t" % "".join([use_dir,outfile])
        self.script += "-l %s \n\n" % "".join([use_dir,logfile])

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


        # Store location of the phylogenetic tree:
        self.sample_data["phylotree"] = self.base_dir + outfile;
        

        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            