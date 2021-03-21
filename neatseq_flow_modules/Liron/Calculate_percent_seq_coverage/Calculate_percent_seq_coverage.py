# -*- coding: UTF-8 -*-
""" 
``Calculate_percent_seq_coverage``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~
    Module to Calculate the percent cover (~alignment) of each of a set of reference sequences within each sample\'s BAM file

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, a sorted bam file type [e.g. mapping result] in:
        ``sample_data[sample]["bam"]``
    * In Project scope a nucleotide file in:
        ``sample_data["project_data"]["fasta.nucl"]``
    * In Sample scope for each Sample, a nucleotide file in:
        ``sample_data[sample]["fasta.nucl"]``    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Calculate_percent_seq_coverage  # Name of the module to use
        base:                                   # Name of the step [or list of names] to run after [must be after a bam file generator step like mapping (must be sorted)]
        script_path:                            # Not needed
        scope:                                  # project or sample the default is project.
        env:                                    # env parameters that needs to be in the PATH for running this module
        use_project_fasta                       # If scope is sample use the project level fasta for all samples
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        redirects:
            -t:                                 # The minimum number of times a given sequence\'s basepair (position) need to be hit by the sample\'s reads, in order to be considered "a covered position"

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"

class Step_Calculate_percent_seq_coverage(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        if "scope" not in list(self.params.keys()):
            self.params["scope"] = 'project'
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
            
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            assert "bam" in list(self.sample_data[sample].keys()), \
                    "In %s:\tThere are no mapping results (bam) for sample %s.\n" % (self.get_step_name(), sample)
            if self.params["scope"] != 'project':
                # Testing for existance of fasta nucleotide file 
                if "use_project_fasta" in list(self.params.keys()):
                    pass
                else:
                    assert "fasta.nucl" in list(self.sample_data[sample].keys()), \
                        "In %s:\tThere is no nucleotide fasta file (fasta.nucl) for sample %s.\n" % (self.get_step_name(), sample)
        if self.params["scope"] == 'project':
            # Testing for existance of fasta nucleotide file 
            assert "fasta.nucl" in list(self.sample_data["project_data"].keys()), \
                "In %s:\tThere is no project level nucleotide fasta file (fasta.nucl) \n" % self.get_step_name()
        else:
            if "use_project_fasta" in list(self.params.keys()):
                assert "fasta.nucl" in list(self.sample_data["project_data"].keys()), \
                    "In %s:\tThere is no project level nucleotide fasta file (fasta.nucl) \n" % self.get_step_name()
        assert "-t" in list(self.params["redir_params"].keys()), \
                "In %s:\tYou must specify a -t redirects option !! \n" % self.get_step_name()
        pass
        
        

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        if "project" in self.params["scope"]:
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()
        
    def build_scripts_byproject(self):
        
        

        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        if (self.params['script_path']==None) and ("multiQueryMultiWGAlign_coverage.py" in os.listdir(self.module_location)):
            self.params['script_path'] = 'python ' + os.path.join(self.module_location,"multiQueryMultiWGAlign_coverage.py")
        
        self.script += "echo '"
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            self.script += "%s\t%%s\n"   % sample % self.sample_data[sample]['bam']
        self.script = self.script.rstrip('\n')
        self.script += "' > %s \n\n" % os.path.join(use_dir,'Samples') 
        
             
        # Define output filename 
        output_filename = "".join([self.sample_data["Title"],'_percent_seq_coverage.tab' ])

        
        self.script += self.get_script_const()
        self.script += "-s %s \\\n\t"   % os.path.join(use_dir,'Samples')
        self.script += "-r %s \\\n\t"   % self.sample_data["project_data"]["fasta.nucl"]
        self.script += "-o %s \n\n"     % os.path.join( use_dir,output_filename)
        
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

        self.sample_data["project_data"]["coverage"] = os.path.join( self.base_dir,output_filename)
        
        self.create_low_level_script()
    
        
        
        
    def build_scripts_bysample(self):
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        if (self.params['script_path']==None) and ("multiQueryMultiWGAlign_coverage.py" in os.listdir(self.module_location)):
            self.params['script_path'] = 'python ' + os.path.join(self.module_location,"multiQueryMultiWGAlign_coverage.py")
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            self.script += "echo '"
            self.script += "%s\t%%s"   % sample % self.sample_data[sample]['bam']
            self.script += "' > %s \n\n" % os.path.join(use_dir,'Samples') 
            
            
            # Define output filename 
            output_filename = "".join([sample,'_', 'percent_seq_coverage.tab'])

            self.script += self.get_script_const()
            
            
            self.script += "-s %s \\\n\t"   % os.path.join(use_dir,'Samples')
            if "use_project_fasta" in list(self.params.keys()):
                self.script += "-r %s \\\n\t"   % self.sample_data["project_data"]["fasta.nucl"]
            else:
                self.script += "-r %s \\\n\t"   % self.sample_data[sample]["fasta.nucl"]
            self.script += "-o %s \n\n"     % os.path.join( use_dir,output_filename)
       
            
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.sample_data[sample]["coverage"] = os.path.join( self.base_dir,output_filename)
            
            self.create_low_level_script()


