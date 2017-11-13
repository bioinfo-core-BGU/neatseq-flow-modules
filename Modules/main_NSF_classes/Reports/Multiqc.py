# -*- coding: UTF-8 -*-
""" 
``Multiqc`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for preparing a MultiQC report for all samples

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* No real requirements. Will give a report with information if one of the base steps produces reports that MultiQC can read, *e.g.* fastqc, bowtie2, samtools etc.
    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* puts report dir in the following slot:

    * ``self.sample_data[<sample>]["Multiqc_report"]``
        
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    firstMultQC:
        module: Multiqc
        base:
        - sam_bwt2_1
        - fqc_trim1
        script_path: /path/to/multiqc

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ewels, P., Magnusson, M., Lundin, S. and Käller, M., 2016. **MultiQC: summarize analysis results for multiple tools and samples in a single report**. *Bioinformatics*, 32(19), pp.3047-3048.    
"""

import os
import sys
import re
from PLC_step import Step,AssertionExcept


__author__ = "Levinl based on Menachem Sklarz"
# Run Multiqc on all base directories

class Step_Multiqc(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = ""
        
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
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
            
            
        #Multiqc main command
        self.script += self.get_script_const()
        # dir_base_list = get_all_bases_dir(self,dir_base_list)
        dir_base_list = " ".join(get_base_dirs_str(self))
        
        self.script += "%s \\\n\t" % dir_base_list
        if "modules" in self.params.keys():
            if self.params["modules"]!="All":
                modules=self.params["modules"].split(",")
                for module in modules:
                    self.script += "-m %s \\\n\t" % module
        elif "modules_exclude" in self.params.keys():
            modules=self.params["modules_exclude"].split(",")
            for module in modules:
                self.script += "-e %s \\\n\t" % module
        self.script += " -s -f  -o %s \n\n" % use_dir
        
        self.sample_data["Multiqc_report"] = self.base_dir
       
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()


        
def get_base_dirs_str(step):
    """ Recursion. Beware
    Returns list of base dirs of the base steps of 'step'
    """

    
    depend_dir_list = []
    
    if step.base_step_list:
        for base_step in step.base_step_list:
            depend_dir_list += ([base_step.base_dir] + get_base_dirs_str(base_step))

    return list(set(depend_dir_list))

    