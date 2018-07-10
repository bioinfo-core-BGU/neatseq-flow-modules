# -*- coding: UTF-8 -*-
""" 
``qiime_filter_otus``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's ``filter_otus_from_otu_table.py``


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A BIOM table in:

    * ``sample_data["biom_table"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting BIOM table in: 

    * ``self.sample_data["biom_table"]``
    
* Puts the BIOM table summary in:

    * ``self.sample_data["biom_table_summary"]``

* Puts the BIOM table in tab-delimited format in:

    * ``self.sample_data["biom_table_tsv"]``

* Puts the unfiltered BIOM table in: 

    * ``self.sample_data["prefilter_biom_table"]``
    

    
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "skip_summary", "", "If passed, will not create the BIOM table summary."
    "skip_tsv", "", "If passed, will not create the tsv version of the BIOM table."
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_filt_otus_1:
        module: qiime_filter_otus
        base: filt_samp_1
        script_path: '{Vars.qiime_path}/filter_otus_from_otu_table.py'
        setenv: {Vars.qiime_env}
        redirects:
            --min_count_fraction: 0.00005
            --min_samples: 10


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


class Step_qiime_filter_otus(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_filter_otus"
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like this is the first qiime step. At the moment, it must come after qiime_prep...\n")
        
        
        
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

        # Defined output files (not final )    
        biom_table = "biom_table.biom"
        filtered_biom_table = "filtered_biom_table.biom"
        biom_table_summary = biom_table + ".summary"
        biom_table_tsv = biom_table + ".txt"
        
        
        ####################################################################################################
        ## 
        ## Step 1: Filtering the biom table
        
        self.script += self.get_script_const()        # Gets the "env", "script_path" and "redir_params" part of the script which is always the same...

        self.script += "-i %s \\\n\t" % self.sample_data["biom_table"]
        self.script += "-o %s \n\n" % "".join([use_dir,biom_table])
        
            
        self.sample_data["prefilter_biom_table"] = self.sample_data["biom_table"]
        self.sample_data["biom_table"] = "".join([self.base_dir,biom_table])
        



        #######################################################################################
        ## 
        ## Step 3: Creating biom table summary
        
        # cmd_text = self.get_script_env_path() 
        if "skip_summary" not in self.params:
            cmd_text = """
    %(script_path)s \\
        -i %(biom)s \\
        -o %(biom_summary)s 
""" % {"biom":"".join([use_dir,biom_table]),\
            "script_path"  : os.sep.join([os.path.split(os.path.normpath(self.params["script_path"]))[0] , \
                                          "biom summarize-table"]),\
            "biom_summary" : "".join([use_dir,biom_table_summary])}
           
            
            self.script += """
# Create summary of biom table for use in rarefaction later

if [ -e %(biom)s ]
then
    %(cmd_text)s
fi
""" % {"biom":"".join([use_dir,biom_table]),\
            "cmd_text":cmd_text}
           
            
            self.sample_data["biom_table_summary"] = "".join([self.base_dir,biom_table_summary])

        ################################################################################################
        ## 
        ## Step 4: Creating biom table in table format
        
        if "skip_tsv" not in self.params:
            cmd_text = """
    %(script_path)s \\
        -i %(biom)s \\
        -o %(biom_tsv)s \\
        --to-tsv \\
        --header-key taxonomy \\
        --output-metadata-id \"Consensus Lineage\"
""" % {"biom":"".join([use_dir,biom_table]),\
            "script_path":os.sep.join([os.path.split(os.path.normpath(self.params["script_path"]))[0] , "biom convert"]),\
            "biom_tsv":"".join([use_dir,biom_table_tsv])}

            self.script += """
# Create biom table in table format

if [ -e %(biom)s ]
then
    %(cmd_text)s
fi

""" % {"biom":"".join([use_dir,biom_table]),\
            "cmd_text":cmd_text}
            
            # Store location of the tsv biom_table:
            self.sample_data["biom_table_tsv"] = "".join([self.base_dir,biom_table_tsv])

        
        
        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

        

        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            