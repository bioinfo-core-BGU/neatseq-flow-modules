# -*- coding: UTF-8 -*-
""" 
``qiime_make_otu_table``
-----------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's ``make_otu_table.py``:

The module creates a BIOM table based on the OTU table and a taxonomy assignment if avaliable (will be available if the ``qiime_assign_taxonomy`` is in the branch).

If chimera checking has been performed, the suspected chimeric sequences will be removed from the BIOM table.

The module also adds code for creating a summary of the BIOM table and a tab-delimited version thereof.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* An OTU table:

    * ``sample_data["otu_table"]``
    
Optional
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A taxonomy assignment of the sequences:

    * ``sample_data["taxonomy"]``
    

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the BIOM table in  

    * ``self.sample_data["project_data"]["biom_table"]``
    
* Puts the BIOM table summary in:

    * ``self.sample_data["project_data"]["biom_table_summary"]``

* Puts the BIOM table in tab-delimited format in:

    * ``self.sample_data["project_data"]["biom_table_tsv"]``

* If a fasta.chimera_removed file exists, will put the unfiltered BIOM table in:

    * ``self.sample_data["project_data"]["unfiltered_biom_table"]``

    
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

    q_mk_otu_1:
        module: qiime_make_otu_table
        base: q_phylo_1
        script_path: '{Vars.qiime_path}/make_otu_table.py'
        setenv: {Vars.qiime_env}
        # skip_summary:
        # skip_tsv:
        redirects:
            --mapping_fp: /path/to/qiime1_mapping.txt


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


class Step_qiime_make_otu_table(Step):
    """ A class that defines a pipeline step name (=instance).
        Most of the class functions are in the super class "Step"
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_make_otu_table"
        
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["project_data"]["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like this is the first qiime step. At the moment, it must come after qiime_prep...\n" )
        
        # If a mapping file was passed, check that the file exists and add it to sample_data
        # This is being done here (and not in step_specific_init(), so that sample_data can be updated as well.
        if "--mapping_fp" in self.params["redir_params"].keys() or "-m" in self.params["redir_params"].keys():
            temp_pf = self.params["redir_params"]["--mapping_fp"] if "--mapping_fp" in self.params["redir_params"].keys() else self.params["redir_params"]["-m"]
            # if not os.path.exists(temp_pf):
                # raise AssertionExcept("The mapping file specified does not exist!!\n")
            if "qiime.mapping" in self.sample_data:
                self.write_warning("Overriding mapping file!")
            self.sample_data["project_data"]["qiime.mapping"] = temp_pf
        
        
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
        ## Step 1: Building OTU table:
        
        self.script += self.get_script_const()        # Gets the "env", "script_path" and "redir_params" part of the script which is always the same...

        self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["otu_table"]
        if "taxonomy" in self.sample_data.keys():
            self.script += "-t %s \\\n\t" % self.sample_data["project_data"]["taxonomy"]
        self.script += "-o %s \n\n" % "".join([use_dir,biom_table])
        
            
        self.sample_data["project_data"]["biom_table"] = "".join([self.base_dir,biom_table])
        
            
            
        ####################################################################################################
        ## 
        ## Step 2: filter out unaligned OTUs AND chimeral sequences (they were removed from the alignment fasta file)
    
        # Do this only if a chimera removal step weas used:
        if "fasta.chimera_removed" in self.sample_data.keys():
            # self.script += self.get_setenv_part()
            # Assuming filter_otus_from_otu_table.py is in the same location as make_otu_table.py used above...
            self.script += "\n\n# Adding code for removal of chimeric sequences from the BIOM table:\n"
            self.script += "%s \\\n\t" % os.sep.join([os.path.split(os.path.normpath(self.params["script_path"]))[0] , "filter_otus_from_otu_table.py"])
            self.script += "-i %s \\\n\t" % (use_dir + biom_table) 
            self.script += "--otu_ids_to_exclude_fp %s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"] 
            self.script += "--negate_ids_to_exclude \\\n\t"
            self.script += "-o %s \n\n" % "".join([use_dir,filtered_biom_table])

            # Store location of unfiltered and active biom tables:
            self.sample_data["project_data"]["unfiltered_biom_table"] = self.sample_data["project_data"]["biom_table"]
            self.sample_data["project_data"]["biom_table"] = "".join([self.base_dir,filtered_biom_table])
            # Change working biom_table to the filtered biom_table:
            biom_table = filtered_biom_table
            biom_table_summary = biom_table + ".summary"
            biom_table_tsv = biom_table + ".txt"



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
           
            
            self.sample_data["project_data"]["biom_table_summary"] = "".join([self.base_dir,biom_table_summary])

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
            self.sample_data["project_data"]["biom_table_tsv"] = "".join([self.base_dir,biom_table_tsv])

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        
        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            