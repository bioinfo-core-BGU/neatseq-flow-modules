# -*- coding: UTF-8 -*-
""" 
``cgMLST_and_MLST_typing``
---------------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad. The MLST typing R script was created by Menachem Sklarz & Michal Gordon 

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~
    A module for a MLST and cgMLST Typing

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Blast results after parsing in: 
        ``self.sample_data[<sample>]["blast.parsed"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Typing results in:
        ``self.sample_data[<sample>]["Typing"]``
    * Merge of typing results in: 
        ``self.sample_data["project_data"]["Typing"]``
    * Files for phyloviz in:
        ``self.sample_data["project_data"]["phyloviz_MetaData"]``
        ``self.sample_data["project_data"]["phyloviz_Alleles"]``
    * Tree file (if --Tree flag is set) in newick format in: 
        ``self.sample_data["project_data"]["newick"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "cut_samples_not_in_metadata",  "", "In the final merge file consider only samples found in the Meta-Data file"
    "sample_cutoff","[0-1]","In the final merge file consider only samples that have at least this fraction of identified alleles"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  The following python packages are required:
    * ``pandas``
*  The following R packages are required:
    * ``magrittr``
    * ``plyr``
    * ``optparse``
    * ``tools``

.. Note:: If using conda environment with R installed the R packages will be automatically installed inside the environment.
        
        
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                   # Name of this step
        module: cgMLST_and_MLST_typing           # Name of the module to use
        base:                                    # Name of the step [or list of names] to run after [must be after steps that generates blast.parsed File_Types] 
        script_path:                             # Leave blank
        metadata:                                # Path to Meta-Data file
        metadata_samples_ID_field:               # Column name in the Meta-Data file of the samples ID
        cut_samples_not_in_metadata:             # In the final merge file consider only samples found in the Meta-Data file
        sample_cutoff:                           # In the final merge file consider only samples that have at least this fraction of identified alleles
        Tree:                                    # Generate newick Tree using hierarchical-clustering [Hamming distance]
        Tree_method:                             # The hierarchical-clustering linkage method [default=complete]
        redirects:
            --scheme:                            # Path to the Typing scheme file [Tab delimited]
            --Type_col_name:                     # Column/s name/s in the scheme file that are not locus names
            --ignore_unidentified_alleles        # Remove columns with unidentified alleles [default=False]

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"

class Step_cgMLST_and_MLST_typing(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ".MLST.type"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        
        assert "--scheme" in self.params["redir_params"] or "-s" in self.params["redir_params"], \
            "In %s:\tYou must redirect a '--scheme' parameter.\n" % self.get_step_name()
        
          
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
            
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
             # Testing for existance of parsed blast data
            assert "blast.parsed" in self.sample_data[sample].keys(), \
                "In %s:\tThere are no parsed blast results for sample %s.\n" % (self.get_step_name(), sample)

        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        # Make a merge file of all results:
        if "Merge_tab_files.py" in os.listdir(self.module_location):
            # Make a dir for the merge file:
            merge_dir = self.make_folder_for_sample("merge")         
            #search for MLST typing files
            look_in=os.sep.join([self.pipe_data["data_dir"].rstrip(os.sep),\
                            self.step , \
                            self.name])
            
            #Running the file merge script
            merge_file=os.sep.join([merge_dir.rstrip(os.sep) , "merged_file.tab"])
            self.script = ""
            self.script += "python %s \\\n\t" % os.path.join(self.module_location,"Merge_tab_files.py")
            self.script += " -D %s \\\n\t" % look_in
            self.script += " -O %s \\\n\t" % merge_file
            self.script += " -R %s$ \n\n" % self.file_tag
            
            self.sample_data["project_data"]["Typing"]=merge_file
            
            if "MLST_parser.py" in os.listdir(self.module_location):                
                # Make a dir for the parsed files:
                pars_dir = self.make_folder_for_sample("Data_for_Phyloviz")
                self.script += "python %s \\\n\t" % os.path.join(self.module_location,"MLST_parser.py")
                #Percentage of identified allele cutoff to consider sample [0.0 - 1.0]
                if "sample_cutoff" in self.params.keys():
                    self.script += " -C %s \\\n\t" % self.params["sample_cutoff"]
                if "metadata" in self.params.keys():
                    self.script += " -M %s \\\n\t" % self.params["metadata"]
                    #samples ID field in the metadata file
                    if "metadata_samples_ID_field" in self.params.keys():
                        self.script += " --S_MetaData %s \\\n\t" % self.params["metadata_samples_ID_field"]
                        #Use only samples with metadata information
                        if "cut_samples_not_in_metadata" in self.params.keys():
                            self.script += " --Cut  \\\n\t" 
                self.script += " -F %s \\\n\t" % merge_file
                if "--Type_col_name" in self.params["redir_params"]:
                    #Non allelic fields in the Merged file
                    self.script += " --Non_allelic '%s,%%s' \\\n\t" % 'Samples,Status,Percentage_of_missing_genes' % self.params["redir_params"]["--Type_col_name"]
                    #Fields in the merge file to move to the metadata file
                    self.script += " --Fields '%s,%%s' \\\n\t" % 'Status,Percentage_of_missing_genes' % self.params["redir_params"]["--Type_col_name"]
                if "Tree" in self.params.keys():
                    self.script += " --Tree  \\\n\t" 
                    if "Tree_method" in self.params.keys():
                        self.script += " --Tree_method %s \\\n\t" % self.params["Tree_method"]
                    self.sample_data["project_data"]["newick"]=os.path.join(pars_dir,"Tree.newick")
 
                self.script += " -O %s \n\n" % pars_dir
                self.sample_data["project_data"]["phyloviz_MetaData"]=os.path.join(pars_dir,"New_MetaData.tab")
                self.sample_data["project_data"]["phyloviz_Alleles"]=os.path.join(pars_dir,"New_Merged_cut.tab")
                
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            if self.params["script_path"]==None:   
                if "MLST_typing.R" in os.listdir(self.module_location):      
                    self.params["script_path"]= "Rscript %s " % os.path.join(self.module_location,"MLST_typing.R")
            
            if self.params["script_path"]!=None:  
                # Define output filename 
                output_filename = "".join([use_dir , sample , self.file_tag])
                
                self.script += self.get_script_const()
                self.script += "--blast %s \\\n\t" % self.sample_data[sample]["blast.parsed"]                
                self.script += "--output %s\n\n" % output_filename
                
                # Store Typing result file:
                self.sample_data[sample]["Typing"] = "".join([sample_dir , sample , self.file_tag])


            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                      
            
            self.create_low_level_script()
                    

                    

def set_global_Sample_data_dir(self,category,info,data):
    if category not in self.sample_data.keys():
        self.sample_data["project_data"][category] = {}
    if self.name not in self.sample_data["project_data"][category].keys():
        self.sample_data["project_data"][category][self.name] = {}
    if info not in self.sample_data["project_data"][category][self.name].keys():
        self.sample_data["project_data"][category][self.name][info] = {}
    self.sample_data["project_data"][category][self.name][info] = data