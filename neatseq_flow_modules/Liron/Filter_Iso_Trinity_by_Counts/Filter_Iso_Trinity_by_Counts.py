# -*- coding: UTF-8 -*-
""" 
``Filter_Iso_Trinity_by_Counts``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~
    A module to preform:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Search for count data in :
        `self.sample_data["project_data"]["isoform.raw_counts"]`
        `self.sample_data["project_data"]["fasta.nucl"]`
    * If 
        
Output:
~~~~~~~~~~~~~

    * Stores Filtered Transcriptome Fasta file in:   
        
        `self.sample_data["project_data"]["fasta.nucl"]`
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  The following R packages are required:
        ``Biostrings``
        ``stringr``
        ``tibble``
        ``tidyverse``
        ``optparse``
        ``ggplot2``
        
        
        
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                      # Name of this step
        module: Filter_Iso_Trinity_by_Counts                              # Name of the used module
        base:                                       # Name of the step [or list of names] to run after with count results.
        script_path:                                # Command for running the a Filter_Iso_Trinity_by_Counts script
                                                    # If this line is empty or missing it will try using the module's associated script
        
        redirects:
            --grouping                              # Path to table of sample groupings. The First column is the samples names (identical to those in --counts header) and a Grouping column ( indicate the name of this column in --grouping_Field  )
            --grouping_Field                        # Name of the Field (Column Name) in the Grouping File (--grouping) to Group the Samples
            --min_count                             # Minimum value of counts table to consider the gene existent in sample.
            --min_sample                            # Minimum samples per group to consider existent in group
            --min_groups                            # Minimum number of groups to consider the gene existent.
            --gene_map                              # Path to gene_map file. If passed, will return representative transcript per gene.
            --plot                                  # Plot an histogram of the counts data over all samples and genes.
            
"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"
# Collect and mereg/append results from all base directories

class Step_Filter_Iso_Trinity_by_Counts(Step):

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        if ("--counts" in self.params["redir_params"]) or ("-c" in self.params["redir_params"] ) :
            raise AssertionExcept("Please do not define  --counts")
        if ("--FASTA" in self.params["redir_params"] ) or ("-f" in self.params["redir_params"] ):
            raise AssertionExcept("Please do not define  --FASTA")
        if ("--output" in self.params["redir_params"] ) or ("-o" in self.params["redir_params"] ) :
            raise AssertionExcept("Please do not define  --output")
        if ("--grouping" not in self.params["redir_params"] ) or ("-g" not in self.params["redir_params"] ) :
            raise AssertionExcept("You Must specify a Grouping File using the --grouping Redirect")
        if ("--grouping_Field" not in self.params["redir_params"] ) or ("-d" not in self.params["redir_params"] ) :
            raise AssertionExcept("You Must specify a Grouping Field using the --grouping_Field Redirect")
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        if 'fasta.nucl' not in list(self.sample_data["project_data"].keys()):
            raise AssertionExcept("No Project level FASTA Nucleotide File was Found")
        if 'isoform.raw_counts' not in list(self.sample_data["project_data"].keys()):
            raise AssertionExcept("No Project level Isoform Raw Counts File was Found")
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        pass
        
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        sample_dir=self.base_dir
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(sample_dir)
        
        if self.params["script_path"]==None: 
            if "filter_trinity_by_counts.R" in os.listdir(self.module_location): 
                self.params["script_path"]= "Rscript  %s "  % os.path.join(self.module_location,"filter_trinity_by_counts.R")
            else:
                raise AssertionExcept("The file %s is not found in the DeSeq2 module directory" % "filter_trinity_by_counts.R" )

        if self.params["script_path"]!=None:    
            # Get constant part of script:
            self.script += self.get_script_const()
            self.script += "--counts %s \\\n\t" % self.sample_data["project_data"]["isoform.raw_counts"]
            self.script += "--FASTA %s \\\n\t"  % self.sample_data["project_data"]["fasta.nucl"]
            self.script += "--output %s \n\n"   % use_dir
            
            
            
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
