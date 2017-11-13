# -*- coding: UTF-8 -*-
""" 
``qiime_divers``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's ``core_diversity_analyses.py``:

The module creates a BIOM table based on the OTU table and a taxonomy assignment if avaliable (will be available if the ``qiime_assign_taxonomy`` is in the branch).

If chimera checking has been performed, the suspected chimeric sequences will be removed from the BIOM table.

The module also adds code for creating a summary of the BIOM table and a tab-delimited version thereof.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A BIOM table:

    * ``sample_data["biom_table"]``
    
Optional
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A phylogenetic tree:

    * ``sample_data["phylotree"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the core diversity directory name in  

    * ``self.sample_data["diversity"]``
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "--mapping_fp", "", "A path to the qiime mapping file (if not set, will use the mapping file passed in ``qiime_prep``."
    "--parameter_fp", "", "A path to a qiime parameter file."

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_divers_1:
        module: qiime_divers
        base: q_filt_otus_1
        script_path: /path/to/QIIME/bin/core_diversity_analyses.py
        qsub_params:
            -pe: shared 20
        sampling_depth: 109897
        redirects:
            --categories: Disease,sex
            --parameter_fp: /path/to/parameter_file


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.




"""



import os
import sys
import re
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_qiime_divers(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_divers"
        
        # Testing the existance of a legitimate parameter file in the pipeline parameter file:
        if not "--parameter_fp" in self.params["redir_params"].keys() and not "-p" in self.params["redir_params"].keys():
            raise AssertionExcept("You must supply a parameter file\n")
        else:
            temp_pf = self.params["redir_params"]["--parameter_fp"] if "--parameter_fp" in self.params["redir_params"].keys() else self.params["redir_params"]["-p"]
            # if not os.path.exists(temp_pf):
                # self.write_warning("Note! The parameter file does not exist\n")
                
                
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like qiime_demult is the first qiime step. At the moment, it must come after qiime_prep...\n" )
        
        if not "sampling_depth" in self.params.keys() and \
           not "--sampling_depth" in self.params["redir_params"] and \
           not "-e" in self.params["redir_params"]:
            raise AssertionExcept("You must supply sampling depth for diversity analysis with redirected parameter --sampling_depth!! The location of the biom table is:\n%s\n" % self.sample_data["qiime"]["biom_table_summary"])
        # This is for backwards compatibility:
        if "sampling_depth" in self.params.keys():
            self.params["redir_params"]["--sampling_depth"] = self.params["sampling_depth"]

        # Testing the existance of a legitimate mapping file in the pipeline parameter file or in self.sample_data:
        # Check if mapping file exists in parameters (overrides mapping from sample_data)
        if "--mapping_fp" in self.params["redir_params"].keys() or "-m" in self.params["redir_params"].keys():
            # Check if mapping file exists in sample_data
            if "qiime.mapping" in self.sample_data.keys():
                self.write_warning("Overriding existing mapping file. Make sure this is OK")
                # mapping_fp = self.sample_data["qiime.mapping"]

            self.sample_data["qiime.mapping"] = self.params["redir_params"]["--mapping_fp"] if "--mapping_fp" in self.params["redir_params"].keys() else self.params["redir_params"]["-m"]
        else:
            if "qiime.mapping" not in self.sample_data.keys():
                raise AssertionExcept("No mapping file exists nor was it passed with -m")

        # # Check if mapping_fp was defined:
        # try:
            # mapping_fp
        # except NameError:
            # raise AssertionExcept("You must supply a mapping file\n")
        # else:
            # pass
            # # Check if mapping_fp exists
            # if not os.path.exists(mapping_fp):
                # self.write_warning("Warning: Mapping file %s does not exist. Is this OK?\n" % mapping_fp)
        
        
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

        
        ####################################################################################################
        ## 
        ## Step 1: Core diversity analyses:
	
	
        # # Defined output files (not final )    
        # biom_table = "biom_table.biom"
        # filtered_biom_table = "filtered_biom_table.biom"
        # biom_table_summary = biom_table + ".summary"
        # biom_table_tsv = biom_table + ".txt"
        
        self.script += self.get_script_const()
        # Add tree if it exists:
        if "phylotree" in self.sample_data.keys():
            self.script += "--tree_fp %s \\\n\t" % self.sample_data["phylotree"]
        # Add mapping file if not passed as redir parameter
        if not "--mapping_fp" in self.params["redir_params"].keys() and not "-m" in self.params["redir_params"].keys():
            self.script += "--mapping_fp %s \\\n\t" % self.sample_data["qiime.mapping"]
        self.script += "-i %(input)s \\\n\t-o %(output)s  \n\n" % {"input" : self.sample_data["biom_table"],\
                                                                                     "output" : use_dir + "core_div/"}
                                                                                    
        self.sample_data["diversity"] = self.base_dir

        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

        

        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            