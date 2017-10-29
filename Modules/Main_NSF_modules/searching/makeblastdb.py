

""" Create a blastdb from a fasta file

        
Requires
~~~~~~~~
        
* fastq files in the following slots:

    * ``sample_data[<sample>]["fasta"]["nucl"|"prot"]``

* Or (if 'projectBLAST' is set)
    * ``sample_data["fasta"]["nucl"|"prot"]``

Output
~~~~~~
    
* A BLAST database in the following slots:

    * ``sample_data[<sample>]["blast"]["blastdb"]["nucl"|"prot"]``
    * ``sample_data[<sample>]["blast"]["blastdb"]["nucl_log"|"prot_log"]``

* Or (if 'projectBLAST' is set):

    * ``sample_data["blast"]["blastdb"]["nucl"|"prot"]``
    * ``sample_data["blast"]["blastdb"]["nucl_log"|"prot_log"]``
    
Parameters that can be set:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide or sample fasta slot should be used"
    "-dbtype", "nucl/prot", "This is a compulsory redirected parameter.Helps the module decide which fasta file to use."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    mkblst1:
        module: makeblastdb
        base: trinity1
        script_path: /path/to/bin/makeblastdb
        redirects:
            -dbtype: nucl
        scope: project

    
"""
    

import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_makeblastdb(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".makeblastdb.out"

        # Checking this once and then applying to each sample:
        if not "-dbtype" in self.params["redir_params"]:
            raise AssertionExcept("You must define a -dbtype parameter\n")
        if not self.params["redir_params"]["-dbtype"] in ["nucl","prot"]:
            raise AssertionExcept("-dbtype must be either nucl or prot\n")
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                self.step_sample_initiation_byproject()
            elif self.params["scope"]=="sample":
                self.step_sample_initiation_bysample()
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
            
        
    def step_sample_initiation_bysample(self):

        # Creating holder for output:
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # Make sure a file exists in the sample equivalent to dbtype:
            try:
                # In version 1.0.2, nucl and prot slots have been renamed to fasta.nucl and fasta.prot
                self.dbtype = self.params["redir_params"]["-dbtype"]
                self.sample_data[sample]["fasta." + self.dbtype]
            except KeyError:
                raise AssertionExcept("No file exists in sample for specified -dbtype (%s)\n" % dbtype, sample)
            # # initialize blast and blastdb slots for sample:
            # if not "blast" in self.sample_data[sample].keys():
                # self.sample_data[sample]["blast"] = dict()
            # if not "blastdb" in self.sample_data[sample]["blast"].keys():
                # self.sample_data[sample]["blast"]["blastdb"] = dict()

    def step_sample_initiation_byproject(self):

        # Make sure a file exists in the sample equivalent to dbtype:
        try:
            self.dbtype = self.params["redir_params"]["-dbtype"]
            self.sample_data["fasta." + self.dbtype]
        except KeyError:
            raise AssertionExcept("No file exists in project for specified -dbtype (%s)\n" % self.dbtype)
            
        # # Creating holder for output:
        # if not "blast" in self.sample_data.keys():
            # self.sample_data["blast"] = dict()
        # if not "blastdb" in self.sample_data["blast"].keys():
            # self.sample_data["blast"]["blastdb"] = dict()
                
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        if self.params["scope"]=="project":
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()

            
            
    def build_scripts_bysample(self):

        
        # Prepare a list to store the qsub names of this steps scripts (will then be put in pipe_data and returned somehow)
        self.qsub_names=[]
        
        # dbtype = "fatsa." + self.params["redir_params"]["-dbtype"]
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
        
            # Check that -dbtype has equivalen file in "fasta" stricture
            # Tested in step_sample_initiation_bysample
            # if not dbtype in self.sample_data[sample]["fasta"]:
                # raise AssertionExcept("No file matching the -dbtype you supplied.\n",sample)
        
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # If "local" is set, will do all IO to local folder and then copy everything to self.base_dir
            if "local" in self.params.keys():
                local_dir = "/local/bioinfo/" + "_".join([self.step,self.name,sample,self.pipe_data["run_code"]]) + os.sep
                self.script += "mkdir -p %s \n\n" % local_dir
                use_dir = local_dir 
            else:
                use_dir = sample_dir

 
 
         
            # Define output filename 
            output_filename = "".join([use_dir , sample , self.file_tag])
            blastdb_title = os.path.basename(output_filename)

            self.script += self.get_script_const()
            self.script += "-out %s \\\n\t" % output_filename
            self.script += "-in %s \\\n\t" % self.sample_data[sample]["fasta." + self.dbtype]
            self.script += "-title %s \\\n\t" % blastdb_title
            self.script += "-logfile %s \n\n" % "%s.log" % output_filename

            
                
            self.sample_data[sample]["blastdb." + self.dbtype] = (sample_dir + os.path.basename(output_filename))
            self.sample_data[sample]["blastdb." + self.dbtype + ".log"] = "%s.log" % (sample_dir + os.path.basename(output_filename))

            self.sample_data[sample]["blastdb"] = self.sample_data[sample]["blastdb." + self.dbtype]
            
            # self.stamp_dir_files(sample_dir)
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                       
            
            self.create_low_level_script()

            
            
            
    def build_scripts_byproject(self):

    
    
        
        # Prepare a list to store the qsub names of this steps scripts (will then be put in pipe_data and returned somehow)
        self.qsub_names=[]
        
        # dbtype = self.params["redir_params"]["-dbtype"]
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
    
        # Check that -dbtype has equivalen file in "fasta" stricture
        # Tested in step_sample_initiation_byproject()
        # if not dbtype in self.sample_data["fasta"]:
            # raise AssertionExcept("No file matching the -dbtype you suppplied.\n")
    
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
        self.script = ""

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
             

     
        # Define output filename 
        output_filename = "".join([use_dir , self.sample_data["Title"] , self.file_tag])
        blastdb_title = os.path.basename(output_filename)

        self.script += self.get_script_const()
        self.script += "-out %s \\\n\t" % output_filename
        self.script += "-in %s \\\n\t" % self.sample_data["fasta." + self.dbtype]
        self.script += "-title %s \\\n\t" % blastdb_title
        self.script += "-logfile %s.log \n\n" % output_filename

        
            
        self.sample_data["blastdb." + self.dbtype] = (self.base_dir + os.path.basename(output_filename))
        self.sample_data["blastdb." + self.dbtype + ".log"] = "%s.log" % (self.base_dir + os.path.basename(output_filename))
        
        self.sample_data["blastdb"] = self.sample_data["blastdb." + self.dbtype]
        # self.stamp_dir_files(self.base_dir)
            
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                    
        
        self.create_low_level_script()

