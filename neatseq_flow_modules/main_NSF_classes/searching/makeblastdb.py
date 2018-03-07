# -*- coding: UTF-8 -*-
""" 
``makeblastdb`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Create a blastdb from a fasta file

        
Requires
~~~~~~~~
        
* fastq files in the following slots:

    * ``sample_data[<sample>]["fasta.nucl"|"fasta.prot"]``

* Or (if 'projectBLAST' is set)
    * ``sample_data["fasta.nucl"|"fasta.prot"]``

Output
~~~~~~
    
* A BLAST database in the following slots:

    * ``sample_data[<sample>]["blastdb"]``
    * ``sample_data[<sample>]["blastdb.nucl"|"blastdb.prot"]``
    * ``sample_data[<sample>]["blastdb.nucl.log"|"blastdb.prot.log"]``

* Or (if 'projectBLAST' is set):

    * ``sample_data["blastdb"]``
    * ``sample_data["blastdb.nucl"|"blastdb.prot"]``
    * ``sample_data["blastdb.nucl.log"|"blastdb.prot.log"]``
    
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

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. **Gapped BLAST and PSI-BLAST: a new generation of protein database search programs**. *Nucleic acids research*, 25(17), pp.3389-3402.
    
"""
    

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_makeblastdb(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "makeblastdb.out"

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
                raise AssertionExcept("No file exists in sample for specified -dbtype (%s)\n" % self.dbtype, sample)

    def step_sample_initiation_byproject(self):

        # Make sure a file exists in the sample equivalent to dbtype:
        try:
            self.dbtype = self.params["redir_params"]["-dbtype"]
            self.sample_data["fasta." + self.dbtype]
        except KeyError:
            raise AssertionExcept("No file exists in project for specified -dbtype (%s)\n" % self.dbtype)
  
  
        
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
        
        
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

         
            # Define output filename 
            output_filename = ".".join([sample, self.name, self.file_tag])
            blastdb_title = os.path.basename(output_filename)

            self.script += self.get_script_const()
            self.script += "-out {dir}{fn} \\\n\t".format(dir= use_dir, fn=output_filename)
            self.script += "-in %s \\\n\t" % self.sample_data[sample]["fasta." + self.dbtype]
            self.script += "-title %s \\\n\t" % blastdb_title
            self.script += "-logfile %s \n\n" % "%s.log" % output_filename

            
                
            self.sample_data[sample]["blastdb." + self.dbtype] = (sample_dir + output_filename)
            self.sample_data[sample]["blastdb." + self.dbtype + ".log"] = "{dir}{fn}.log".format(dir= sample_dir, fn=output_filename)

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
        output_filename = ".".join([self.sample_data["Title"], self.name, self.file_tag])
        blastdb_title = os.path.basename(output_filename)

        self.script += self.get_script_const()
        self.script += "-out {dir}{fn} \\\n\t".format(dir= use_dir, fn=output_filename)
        self.script += "-in %s \\\n\t" % self.sample_data["fasta." + self.dbtype]
        self.script += "-title %s \\\n\t" % blastdb_title
        self.script += "-logfile %s.log \n\n" % output_filename

        
            
        self.sample_data["blastdb." + self.dbtype] = (self.base_dir + output_filename)
        self.sample_data["blastdb." + self.dbtype + ".log"] = "{dir}{fn}.log".format(dir= self.base_dir, fn=output_filename)
        
        self.sample_data["blastdb"] = self.sample_data["blastdb." + self.dbtype]
        # self.stamp_dir_files(self.base_dir)
            
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                    
        
        self.create_low_level_script()

