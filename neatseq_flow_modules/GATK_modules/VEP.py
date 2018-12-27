#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

""" A class defining a pipeline.

This class takes input files: samples and parameters, and creates a qsub pipeline, including dependencies
Actual work is done by calling other class types: PLCStep and PLCName
"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_VEP(Step):
    """ A class that defines a pipeline step name (=instance).
    """


    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
            
            
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Prepare a list to store the qsub names of this steps scripts (will then be put in pipe_data and returned somehow)
        self.qsub_names=[]
        
       
        # Each iteration must define the following class variables:
        # spec_qsub_name
        # spec_script_name
        # script
######################################################## SNP

        use_dir = self.local_start(self.base_dir)
        for chr in self.params["chrom_list"].split(','):
            chr = chr.strip()
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],chr])

            input_file = self.sample_data[chr]["vcf"]

            output_file = use_dir + self.sample_data["Title"] + "_VEP_" + chr + ".vcf"
            self.script = ""
            self.script += self.get_script_const()
            self.script += """-i %(input_file)s \\
            -o %(output_file)s\n
            """ % { 
                    "input_file" : input_file,
                    "output_file" : output_file
            }      
        

            self.sample_data[chr]["vcf"]= use_dir + self.sample_data["Title"] + "_VEP_" + chr + ".vcf"
            self.stamp_file(self.sample_data[chr]["vcf"])
            
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
                        

