# -*- coding: UTF-8 -*-
""" 
``Qualimap``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~
    Runs Qualimap in RNASeq or BAM qc modes, in RNASeq mode will perform count analysis as well 

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, a bam file in:
        ``self.sample_data[sample]["bam"]``
    * For RNASeq mode a project level GTF file in (or can be passed using the -gtf redirects argument):
        ``self.sample_data["project_data"]["GTF"]``
        

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, a directory with Qualimap results in:
        ``self.sample_data[sample]["qualimap"]``
    * In RNASeq mode for each Sample, a count data in:
        ``self.sample_data[sample]["HTSeq.counts"]``
        

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "mode",  "rnaseq,bamqc", "Qualimap mode of action"
    "count_threshold",  "int", "Threshold for the number of counts"
    
Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Qualimap                        # Name of the module to use
        base:                                   # Name of the step [or list of names] to run after [must be after a bam generating step]
        script_path:                            # Command for running Qualimap 
        setenv:                                 # env parameters that needs to be in the PATH for running this module
        mode:                                   # Qualimap mode of action
        gff2infofile:                           # For rnaseq mode, python script to convert GFF to info file, 
                                                # by default will use the Qualimap's script (from https://bitbucket.org/kokonech/qualimap/raw/deb106fabc8c5c58319beab4f659b1c3a0bffc60/util/createQualimapInfoFile.py )
        genome:                                 # For rnaseq mode, a reference genome file (FASTA)
        count_threshold:                        # Threshold for the number of counts 
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        redirects:
            -gtf:                               # For rnaseq mode, annotations file in Ensembl GTF format.
            --java-mem-size=20G:                # Parameters for setting the amount of memory available for java
            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2015). “Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data.” Bioinformatics, btv566.‏
"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"

class Step_Qualimap(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ""
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        assert "mode" in list(self.params.keys()), \
            "In %s:\tA mode of action must be specify!!.\n" % (self.get_step_name())
            
        if (self.params["mode"]=='rnaseq') and ('-gtf' not in list(self.params["redir_params"].keys())):
            assert "GTF" in list(self.sample_data["project_data"].keys()), \
                "In %s:\tA GTF file must be be specified using the -gtf redirects argument or in a GTF project level slot.\n" % (self.get_step_name())
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # Testing for existance of binning data
            assert "bam" in list(self.sample_data[sample].keys()), \
                "In %s:\tThere is no bam file for sample %s.\n" % (self.get_step_name(), sample)
        pass
        

    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        #initiating new script 
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
       
        if self.params["mode"]=='rnaseq':
            
            Count_Analysis = self.make_folder_for_sample("Count_Analysis")
            use_Count_Analysis = self.local_start(Count_Analysis)
            self.script += "cd %s \n\n" % use_Count_Analysis
            self.script = "echo -e '#Sample\\tCondition\\tPath\\tColumn\\n"
            for sample in self.sample_data["samples"]:
                self.script += " %s\\t1\\t%%s\\t2\\n" % sample % self.sample_data[sample]["HTSeq.counts"]
            
            self.sample_data["project_data"]["count_info"]=os.path.join(Count_Analysis,'count_info')
            self.script += "' > %s\n\n" % self.sample_data["project_data"]["count_info"]
            self.sample_data["project_data"]["count_info"]=os.path.join(use_Count_Analysis,'count_info')
            
            if (len(self.sample_data["samples"])>10):
                self.script += "script_dir=$(which %s) \n\n" % self.params["original_script_path"]
                self.script += "script_dir=${script_dir%%/*}/scripts/%s \n\n" % 'countsQC.r'
                self.script += "if [ -f  '$script_dir' ]; then\n"
                self.script += "sed -i 's/image.width <- 3\*480/image.width <- 10\*480/g' $script_dir ; \n" 
                self.script += "sed -i 's/image.height <- 3\*480/image.height <- 10\*480/g' $script_dir ; \n" 
                self.script += "fi\n\n"
                
            if "gff2infofile" not in list(self.params.keys()):
                self.params["gff2infofile"]=None
            if self.params["gff2infofile"]==None:
                if "createQualimapInfoFile.py" in os.listdir(self.module_location):
                    self.params["gff2infofile"]=os.path.join(self.module_location,"createQualimapInfoFile.py")
            
            self.sample_data["project_data"]["infofile"]=None
            if (self.params["gff2infofile"]!=None) and ("genome" in list(self.params.keys())):
                if self.params["genome"]!=None:
                    self.sample_data["project_data"]["infofile"]=os.path.join(use_Count_Analysis,'annotations.txt')
                    self.script += "python %s \\\n\t" % self.params["gff2infofile"]
                    if '-gtf' in list(self.params["redir_params"].keys()):
                        self.script += "-g %s \\\n\t"  % self.params["redir_params"]["-gtf"]
                    else:
                        self.script += "-g %s \\\n\t"  % self.sample_data["project_data"]["GTF"]
                    self.script += "-f %s \\\n\t"  % self.params["genome"]
                    self.script += "-o %s \n\n" %  self.sample_data["project_data"]["infofile"]
                    
            self.script += "%s counts \\\n\t" % self.params["original_script_path"]
            self.script += "-d %s \\\n\t"     % self.sample_data["project_data"]["count_info"]
            if "count_threshold" in list(self.params.keys()):
                self.script += "-k %s \\\n\t"     % self.params["count_threshold"]
            if self.sample_data["project_data"]["infofile"]!=None:
                self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["infofile"]
                self.sample_data["project_data"]["infofile"]=os.path.join(Count_Analysis,'annotations.txt')
            for mem in list(self.params["redir_params"].keys()):
                if mem.startswith('--java-mem-size'):
                    if self.params["redir_params"][mem]!=None:
                        self.script += "%s%%s%%%%s \\\n\t" % mem % self.params["arg_separator"] % self.params["redir_params"][mem]
                    else:
                        self.script += "%s \\\n\t" % mem
            self.script += "-outformat PDF -outdir  %s \n\n"  % use_Count_Analysis
            
            if (len(self.sample_data["samples"])>10):
                self.script += "script_dir=$(which %s) \n\n" % self.params["original_script_path"]
                self.script += "script_dir=${script_dir%%/*}/scripts/%s \n\n" % 'countsQC.r'
                self.script += "sed -i 's/image.width <- 10\*480/image.width <- 3\*480/g' $script_dir \n\n" 
                self.script += "sed -i 's/image.height <- 10\*480/image.height <- 3\*480/g' $script_dir \n\n" 
            
            
            self.local_finish(use_Count_Analysis,Count_Analysis)       # Sees to copying local files to final destination (and other stuff)
           
            
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        self.build_scripts_bysample()
        
        
    def build_scripts_bysample(self):
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
            
        self.params["original_script_path"] = self.params["script_path"]
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            self.params["script_path"]= self.params["script_path"] + " " + self.params["mode"]
            self.script += self.get_script_const()
            if self.params["mode"]=='rnaseq':
                if '-gtf' not in list(self.params["redir_params"].keys()):
                    self.script += "-gtf %s \\\n\t"  % self.sample_data["project_data"]["GTF"]
                self.sample_data[sample]["HTSeq.counts"] = os.path.join(use_dir,sample+'.counts')
                self.script += "-oc %s \\\n\t" % self.sample_data[sample]["HTSeq.counts"]
                self.sample_data[sample]["HTSeq.counts"] = os.path.join(sample_dir,sample+'.counts')
                if ("fastq.F" in list(self.sample_data[sample].keys())) and ("fastq.R" in list(self.sample_data[sample].keys())):
                    if ('-pe' not in list(self.params["redir_params"].keys())) or ('--paired' not in list(self.params["redir_params"].keys())):
                        self.script += "-pe  \\\n\t"
            self.script += "-bam %s \\\n\t" %  self.sample_data[sample]["bam"]
            self.script += "-outdir %s \\\n\t" %  use_dir
            
            self.script +="\n\n"
            self.create_low_level_script()
