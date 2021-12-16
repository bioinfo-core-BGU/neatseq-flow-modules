# -*- coding: UTF-8 -*-
""" 
``ProteinOrtho``
-----------------------

:Authors: Liron Levin
:Affiliation:  Bioinformatics core facility
:Organization: Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~
    A module for running ProteinOrtho

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * proteins file in 
        ``self.sample_data[sample]["fasta.prot"]``


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * puts output cluster tab file in:
        ``self.sample_data["project_data"]['matrix']``
        
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "Nsplit",  "None", "Number of Splits to partition the data and parallel the analysis"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``ProteinOrtho v``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                                   # Name of this step
        module: ProteinOrtho                                     # Name of the module used
        base:                                                    # Name of the step [or list of names] to run after [must be after a fasta.prot file generator step ]
        script_path:                                             # Command for running the ProteinOrtho script 
        qsub_params:
            -pe:                                                 # Number of CPUs to reserve for this analysis
        Nsplit:                                                  # Number of Splits to partition the data and parallel the analysis
        redirects:
            -project:                                            # Name Of the Project
            -cpus:                                               # Number of CPUs to use in this analysis for each split


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Li, Bo, and Colin N. Dewey. "RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome." BMC bioinformatics 12.1 (2011): 323.‏

"""


import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Levin Levin"
__version__= ""


class Step_ProteinOrtho(Step):
 
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ""
        self.params['arg_separator'] = "="
        if '-project' not in list(self.params["redir_params"].keys()):
            self.params["redir_params"]["-project"] = "ProteinOrtho"
        assert "Nsplit"  in list(self.params.keys()) , \
            "you should  provide a Number of Splits to parallel the analysis [ step %s ]\n" % self.get_step_name()
        assert not ("-jobs" in list(self.params["redir_params"].keys())) , \
            "you can't use '-jobs' option, this will be set by the module. [ step %s ]\n" % self.get_step_name()
        assert not ("-step" in list(self.params["redir_params"].keys())) , \
            "you can't use '-step' option, this will be set by the module. [ step %s ]\n" % self.get_step_name()

        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        for sample in self.sample_data["samples"]:    
            assert "fasta.prot" in self.sample_data[sample].keys(), "Sample %s does not have fasta.prot files in step %s\n if you have bam files use --bam\n" % (sample, self.name)
        pass
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        # Step 1
        self.ProteinOrtho_dir = self.make_folder_for_sample("ProteinOrtho")
        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.ProteinOrtho_dir)
        #initiating new script 
        self.script = "cd %s\n\n" % use_dir
        
        # Get constant part of script:
        self.script += self.get_script_const()
        self.script +="-step=1 \\\n\t"
        for sample in self.sample_data["samples"]:
            self.script += "%s " % self.sample_data[sample]["fasta.prot"]
        
        self.script += "\n\n"
        self.local_finish(use_dir,self.ProteinOrtho_dir) 
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
                # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.ProteinOrtho_dir)
        #initiating new script 
        self.script = "cd %s\n\n" % use_dir
        
        # Get constant part of script:
        self.script += self.get_script_const()
        self.script +="-step=3 \\\n\t"
        for sample in self.sample_data["samples"]:
            self.script += "%s " % self.sample_data[sample]["fasta.prot"]
        
        self.script += "\n\n"
        self.local_finish(use_dir,self.ProteinOrtho_dir) 
        pass
        
                    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        if isinstance(self.params["Nsplit"], str):
            if self.params["Nsplit"].isdigit():
                self.params["Nsplit"] = int(self.params["Nsplit"])
        
        
        if isinstance(self.params["Nsplit"], int):
            for split in range(self.params["Nsplit"]):

                use_dir = self.ProteinOrtho_dir
                self.spec_script_name = self.set_spec_script_name("Step_"+str(split+1))
                #initiating new script 
                self.script = "cd %s\n\n" % use_dir
                
                # Get constant part of script:
                self.script += self.get_script_const()
                self.script +="-step=2 \\\n\t"
                self.script +="-jobs=%s/%%s \\\n\t" % str(split+1) \
                                                    % int(self.params["Nsplit"])
                for sample in self.sample_data["samples"]:
                    self.script += "%s " % self.sample_data[sample]["fasta.prot"]
                self.script += "\n\n"
                self.create_low_level_script()
        self.sample_data["project_data"]['matrix'] = os.sep.join([use_dir.rstrip(os.sep),self.params["redir_params"]['-project']+".proteinortho.tsv"]) 
            
def set_Sample_data_dir(self,category,info,data):
    if category not in list(self.keys()):
        self[category] = {}
    if info not in list(self[category].keys()):
        self[category][info] = {}
    self[category][info] = data 