# -*- coding: UTF-8 -*-
""" 
``Gassst``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~
    A module for executing Gassst on a nucleotide fasta file.
    The search can be either on a sample fasta or on a project-wide fasta.
    It can use the fasta as a database or as a query.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * fasta files in the following slot for sample-wise Gassst:
        * ``sample_data[<sample>]["fasta.nucl"]``
    * or fasta files in the following slots for project-wise Gassst:
        * ``sample_data["fasta.nucl"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * puts Gassst output files in the following slots for sample-wise Gassst:
        * ``sample_data[<sample>]["blast"]``
        * ``sample_data[<sample>]["blast.nucl"]``
    * puts fasta output files in the following slots for project-wise Gassst:
        * ``sample_data["blast"]``
        * ``sample_data["blast.nucl"]``
        
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope",  "project/sample", "Set if project-wide fasta.nucl file type should be used [project] the default is sample-wide fasta.nucl file type"
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``Gassst v1.28``
    *  The following python packages are required:
        ``pandas``
    * Only -d [database] or -i [query] not both
    * The Gassst module will generate blast like output with fields:
        ```"qseqid sallseqid qlen slen qstart qend sstart send length evalue sseq"``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                         # Name of this step
        module: Gassst                 # Name of the module to use
        base:                          # Name of the step [or list of names] to run after [mast be after a fasta generating step]
        script_path:                   # Command for running the Gassst script
                                       # The Gassst module will generate blast like output with fields:
                                       # "qseqid sallseqid qlen slen qstart qend sstart send length evalue sseq"
        scope:                         # Set if project-wide fasta.nucl file type should be used [project] the default is sample-wide fasta.nucl file type
        qsub_params:
            -pe:                       # Number of CPUs to reserve for this analysis
        redirects:
            -h:                        # Max hits per query, for downstream best hit will be chosen!
            -i:                        # Only -d [database] or -i [query] not both
            -l:                        # Complexity_filter off
            -d:                        # Only -d [database] or -i [query] not both
            -n:                        # Number of CPUs running Gassst
            -p:                        # Minimum percentage of identity. Must be in the interval [0 100]

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Rizk, Guillaume, and Dominique Lavenier. "GASSST: global alignment short sequence search tool." Bioinformatics 26.20 (2010): 2534-2540.‏
"""


import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"

class Step_Gassst(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ".Gassst.out"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        # Gassst can only work with nucleotide fasta
        self.params["fasta2use"] = "nucl"
        # Check that either db or query (not both) are set in redir_params:
        assert len(set(["-d","-i"]) & set(self.params["redir_params"].keys())) == 1, "In %s:\tYou must supply either 'db' or 'query'\n" % self.get_step_name()       
        # Check that the -p argument is supplied
        assert "-p" in self.params["redir_params"].keys(), "In %s:\tYou must supply -p argument for Gassst \n" % self.get_step_name()       

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" in self.params.keys():
            if self.params["scope"]=="project":
                self.step_sample_initiation_byproject()
            else:
                self.step_sample_initiation_bysample()
        else:
            self.step_sample_initiation_bysample()
        
        
        
    def step_sample_initiation_bysample(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for sample-level BLAST
        """
        
            
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            if not "blast" in self.sample_data[sample].keys():
                self.sample_data[sample]["blast"] = dict()
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            if not "blast.nucl" in self.sample_data[sample].keys():
                self.sample_data[sample]["blast.nucl"] = dict()

            # Decide on locations of db and query
            if "-i" in self.params["redir_params"].keys():
                assert "fasta.nucl" in self.sample_data[sample].keys(), "In %s:\tFor sample-as-DB , you need to have a fasta in the sample (sample %s).\nIf the query is a project fasta, set parameter 'scope' to 'project' \n" % (self.get_step_name(), sample)
                
            # Decide which fasta to use in Gassst:
            # "fasta" is not defined for the sample:
            assert "fasta.nucl" in self.sample_data[sample].keys(), "In %s:\tNo 'fasta.nucl' defined for sample %s.\nIf the query is a project fasta, use parameter 'scop: project'\n" % (self.get_step_name(),sample)       
        
        pass

    def step_sample_initiation_byproject(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for project-level BLAST
        """
        
            
        if not "blast" in self.sample_data.keys():
            self.sample_data["project_data"]["blast"] = dict()
        
        if not "blast.nucl" in self.sample_data.keys():
            self.sample_data["project_data"]["blast.nucl"] = dict()

        assert "fasta.nucl" in self.sample_data.keys(), "In %s:\tYou need a 'fasta.nucl' file defined to run Gassst.\nIf the 'fasta.nucl' files are per sample, use 'scope: sample' parameter.\n" % (self.get_step_name())
            
            
        # Decide on locations of db and query
        if "-i" in self.params["redir_params"].keys():
            assert "fasta.nucl" in self.sample_data.keys(), "In %s:\tFor sample-as-DB , you need to have a fasta.nucl in the sample  .\n" % (self.get_step_name())
            
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """

        if "scope" in self.params.keys():
            if self.params["scope"]=="project":
                pass
            else:
                self.make_sample_file_index()   # see definition below
        else:
            self.make_sample_file_index()   # see definition below

        
        

    
    def build_scripts(self):
        """ This is the actual script building function
            
        """
        if "scope" in self.params.keys():
            if self.params["scope"]=="project":
                self.build_scripts_byproject()
            else:
                self.build_scripts_bysample()
        else:
            self.build_scripts_bysample()
        

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
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
                
                
            # Define output filename 
            output_filename = "".join([use_dir , sample , self.file_tag])
            if "Gassst2blast.py" in os.listdir(self.module_location):
                    self.params["redir_params"]["-m"]='0'
            self.script += self.get_script_const()
            # Define query and db files:
            # If db is defined by user, set the query to the correct 'fasta'
            if "-d" in self.params["redir_params"].keys():
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["fasta.nucl"]
            # If db is not defined by user, set the db to the correct blastdb, with 'fasta'
            # query must be set by user. assertion is made in step_specific_init()
            else:
                self.script += "-d %s \\\n\t" % self.sample_data[sample]["fasta.nucl"]
                
            self.script += "-o %s\n\n" % output_filename
            
            if "Gassst2blast.py" in os.listdir(self.module_location):
                self.script += "python %s \\\n\t" % os.path.join(self.module_location,"Gassst2blast.py")
                self.script += "-i %s \\\n\t" % output_filename
                self.script += "-o %s \\\n\t" % output_filename
            # Store BLAST result file:
            self.sample_data[sample]["blast"] = (sample_dir + os.path.basename(output_filename))
            self.stamp_file(self.sample_data[sample]["blast"])
            
            self.sample_data[sample]["blast.nucl"] = (sample_dir + os.path.basename(output_filename))
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            
            self.create_low_level_script()
                    
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """


        
        
        # Each iteration must define the following class variables:
        # spec_script_name
        # script
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
                
                
        # Define output filename 
        output_filename = "".join([use_dir , self.sample_data["Title"] , self.file_tag])
        if "Gassst2blast.py" in os.listdir(self.module_location):
            self.params["redir_params"]["-m"]='0'
        self.script += self.get_script_const()
        # Define query and db files:
        # If db is defined by user, set the query to the correct 'fasta'
        if "-d" in self.params["redir_params"].keys():
            self.script += "-i %s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
        # If -d is not defined by user, set the -d to the correct fasta, with 'fasta2use'
        # -i must be set by user. assertion is made in step_specific_init()
        else:
            self.script += "-d %s \\\n\t" % self.sample_data["project_data"]["fasta.nucl"]
                
        self.script += "-o %s\n\n" % output_filename
        
        if "Gassst2blast.py" in os.listdir(self.module_location):
            self.script += "python %s \\\n\t" % os.path.join(self.module_location,"Gassst2blast.py")
            self.script += "-i %s \\\n\t" % output_filename
            self.script += "-o %s \\\n\t" % output_filename
        # Store BLAST result file:
        self.sample_data["project_data"]["blast"] = (self.base_dir + os.path.basename(output_filename))
        self.stamp_file(self.sample_data["project_data"]["blast"])
        self.sample_data["project_data"]["blast.nucl"] = (self.base_dir + os.path.basename(output_filename))

        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
                

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names.
            This can be used by scripts called by create_spec_wrapping_up_script() to summarize the BLAST outputs.
        """
        
        with open(self.base_dir + "Gassst_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tGassst_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["blast"]))
                
        self.sample_data["project_data"]["BLAST_files_index"] = self.base_dir + "Gassst_files_index.txt"
        
  
        
