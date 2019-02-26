# -*- coding: UTF-8 -*-

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

class Step_Picard_CollectAlignmentSummaryMatrics(Step):
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
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample])
            self.script = ""
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            my_string = """
                cd %(sample_dir)s
                echo '\\n---------- CollectAlignmentSummaryMetrics -------------\\n'
                %(picard_path)s CollectAlignmentSummaryMetrics \\
                R=%(genome_reference)s \\
                INPUT=%(BAM)s \\
                OUTPUT=%(output)s \\

            """ % { "sample_dir" : sample_dir,
                    "picard_path" : self.params["script_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "BAM" : self.sample_data[sample]["bam"],
                    "output" : sample_dir + sample + "_CollectAlignmentSummaryMetrics"

            }
            
            
            
            self.script += my_string
            #self.get_script_env_path()
            
            
            
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                    
