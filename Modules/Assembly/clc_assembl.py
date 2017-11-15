# -*- coding: UTF-8 -*-
""" 
``clc_assembl``
--------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for assembling reads using CLC assembler.

 
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* fastq files in at least one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    
Output:
~~~~~~~~~~~~~


* puts fasta output files in the following slots:
    
    * if ``scope`` set to ``sample``:
    
        * ``sample_data[<sample>]["fasta.nucl"]``
        * ``sample_data[<sample>]["clc_assembl.contigs"]``
        
        * Also, sets ``sample_data[<sample>]["assembler"] = "clc_assembl"``
    
    * if ``scope`` set to ``project``:
    
        * ``sample_data["fasta.nucl"]``
        * ``sample_data["clc_assembl.contigs"]``

        * Also, sets ``sample_data[<sample>]["assembler"] = "clc_assembl"``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set to ``project`` to assembl all project reads into one assembly."
    "p", "e.g. 'fb ss 180 250'", "Sets the ``-p`` parameter passed to CLC for paired-end reads. Required only if the project includes paired end reads."
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    clc1:
        module: clc_assembl
        base: trim1
        script_path: /path/to/clc_assembler
        qsub_params:
            -pe:    shared 30
            node:   sge37
        scope:      sample
        p:          fb ss 180 250 
        redirects:
            --cpus: 30
            

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_clc_assembl(Step):

  

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".clc_assembl.fna"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)

        
        if "scope" in self.params:
          
            if self.params["scope"] !="project" and self.params["scope"] != "sample":
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
        
            
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
    
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        if self.params["scope"] == "project":

        
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""

            # Make a dir for the current sample:
            output_file = os.sep.join([self.base_dir, self.sample_data["Title"] + self.file_tag])
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            self.script += self.get_script_const()
            self.script += "-o %s \\\n\t" % output_file

            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                if "PE" in self.sample_data[sample]["type"] and "SE" in self.sample_data[sample]["type"]:
                    print >> sys.stdout, "CLC assembler not defined for PE-SE mixes. Using PE file only..."
                if "PE" in self.sample_data[sample]["type"]:
                    try:
                        self.script += "-p %s \\\n\t-q \\\n\t" % self.params["p"]
                    except KeyError:
                        raise AssertionExcept("With paired end reads, you must specify a 'p' parameter containing information to pass with '-p' to clc_assembler. See the clc manual.")
                    self.script += "-i %s %s \\\n\t" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
                elif "SE" in self.sample_data[sample]["type"]:
                    self.script += "-p no \\\n\t-q \\\n\t%s \\\n\t" % self.sample_data[sample]["fastq.S"]
                else:       # Mixed!!
                    pass
                
            # Remove trailing '\\\n\t' from last iteration and add some newlines for clarity
            self.script = self.script.rstrip("\\\n\t") + "\n\n"

            # Store results to fasta and assembly slots:
            self.sample_data["fasta.nucl"] = output_file
            self.sample_data[self.step + ".contigs"] = output_file
            self.sample_data["assembler"] = self.get_step_step()

            self.stamp_file(self.sample_data[self.step + ".contigs"])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                        
        
        else:
        
            # Each iteration must define the following class variables:
                # spec_script_name
                # script
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)
                output_file = os.sep.join([sample_dir, sample + self.file_tag])
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                self.script += self.get_script_const()
                self.script += "-o %s \\\n\t" % output_file

                if "mixed" in self.sample_data[sample]["type"]:
                    print >> sys.stdout, "CLC assembler not defined for PE-SE mixes. Using PE file only..."
                if "PE" in self.sample_data[sample]["type"] or "mixed" in self.sample_data[sample]["type"]:
                    self.script += "-p %s \\\n\t-q \\\n\t" % self.params["p"]
                    self.script += "-i %s %s\n\n" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
                elif "SE" in self.sample_data[sample]["type"]:
                    self.script += "-p no \\\n\t-q \\\n\t%s\n\n" % self.sample_data[sample]["fastq.S"]
                else:       # Mixed!!
                    pass
                    

                # Store results to fasta and assembly slots:
                self.sample_data[sample]["fasta.nucl"] = output_file
                self.sample_data[sample][self.step + ".contigs"] = output_file
                self.sample_data[sample]["assembler"] = self.get_step_step()

                self.stamp_file(self.sample_data[sample][self.step + ".contigs"])

                    
                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                
                self.create_low_level_script()
                        

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by downstream collection scripts
        """
        
        pass