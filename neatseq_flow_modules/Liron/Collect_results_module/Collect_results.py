#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python
# -*- coding: UTF-8 -*-
""" 
Module ``Collect_results``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

SHORT DESCRIPTION
    A module to Collect and merge/append results from all base steps directories:
    This module will search for each base step for all the results files with a common name pattern [Regular expression].
    The search will be done within the base step result directories.
    The sample name could be inferred for each result file base on the parent directory name and added to the merged file [as new column named "Samples"].
    All the results files will be append [by default] or merged by a common column name.
    The merge files can then be convert individually to pivot table file

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Tab delimited files with common name pattern found within the base step data directories:
    * ``For example files ending with .out``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Generate merged tab delimited files:
    * ``Will generate file for each of the base steps with the file ending with .merg``
    * Can also generate Excel file with sheet for each base step 


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "PARAMETER NAME",  "POSSIBLE VALUES", "DESCRIPTION"
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    Step_Name:                            # Name of this step
        module: Collect_results           # Name of the used module
        base:                             # Name of the step [or list of names] to run after and collect results from [must be after a merge step]
        script_path:                      # Command for running the a merging script
                                          # If this line is empty or missing it will try using the module's associated script
        redirects:
            -R:                           # Regular expression to find result files
            --Merge_by:                   # Merge files by common column
            --header:                     # Don't use a header row, use integers instead [0,1,2,3...], easy to use with --pivot option
            --Excel:                      # Collect all results to excel file split by sheets
            --add_samples_names:          # Infer and add samples names from file parent directory to "Samples" column
            --pivot:                      # Convert to pivot table by [index columns values]
                                          # If with the options: -add_samples_names and --header  it is possible to use: '''Samples'' '5' '0''
            --MetaData:                   # Use external MetaData file as the base for merging
            --split_by:                   # Split the data in the columns [index <columns> values] before pivot
            --sep:                        # Columns separator for input file
            -T:                           # Write Transpose output
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
# Collect and mereg/append results from all base directories

class Step_Collect_results(Step):

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".merg"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        pass
        
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
        self.spec_script_name = "_".join([self.step,self.name])
        self.script = ""
        sample_dir=self.base_dir
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(sample_dir)
        
        if self.params["script_path"]==None: 
            if "Merge_tab_files.py" in os.listdir(self.module_location): 
                self.params["script_path"]= "python  %s "  % os.path.join(self.module_location,"Merge_tab_files.py")
            else:
                raise AssertionExcept("The file %s is not found in the Collect_results module directory" % "Merge_tab_files.py" )
        
        if self.params["script_path"]!=None:    
            #Collect_results main command
            for base in self.get_base_step_list():
                if "env" in self.params.keys():
                    self.script +="env %s  \\\n" % self.params["env"]
                self.script += "%s \\\n\t"  % self.params["script_path"]
                for par in self.params["redir_params"]:
                    if par not in ["-D","-O"]:
                        if self.params["redir_params"][par]!=None:
                            self.script += "%s %%s \\\n\t"  % par \
                                                            % self.params["redir_params"][par]
                        else:
                            self.script += "%s  \\\n\t"  % par 
                self.script += "-D %s \\\n\t" % base.base_dir
                if "--Excel" in self.params["redir_params"]:
                    self.script += "-O %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),"Collected_results.xlsx"])
                else:
                    self.script += "-O %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),base.name+self.file_tag])
        
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
