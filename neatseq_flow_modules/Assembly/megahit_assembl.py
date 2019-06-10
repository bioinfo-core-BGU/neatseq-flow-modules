# -*- coding: UTF-8 -*-
""" 
``megahit_assembl``
--------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for assembling reads using MEGAHIT assembler.

 
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
        * ``sample_data[<sample>]["megahit_assembl.contigs"]``
        
        * Also, sets ``sample_data[<sample>]["assembler"] = "megahit_assembl"``
    
    * if ``scope`` set to ``project``:
    
        * ``sample_data["fasta.nucl"]``
        * ``sample_data["megahit_assembl.contigs"]``

        * Also, sets ``sample_data[<sample>]["assembler"] = "megahit_assembl"``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set to ``project`` to assembl all project reads into one assembly."
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    megahit1:
        module: megahit_assembl
        base: trim1
        script_path: /path/to/megahit
        qsub_params:
            -pe: shared 30
            node: sge37
        scope: project
        redirects:
            --continue: 
            --num-cpu-threads: 30

            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Li, D., Liu, C.M., Luo, R., Sadakane, K. and Lam, T.W., 2015. **MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph**. *Bioinformatics*, 31(10), pp.1674-1676.


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_megahit_assembl(Step):



    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".megahit.out"

        if "scope" not in list(self.params.keys()):
            self.params["scope"] = "sample"
            self.write_warning("You did not specify 'scope'. Setting to 'sample'\n")
            

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)
            if len({"fastq.F", "fastq.R"} & set(self.sample_data[sample].keys())) ==1:
                raise AssertionExcept("Sample has only forward or reverse reads. It must have either pairs or single reads\n", sample)

        
        if "scope" in self.params:
          
            if self.params["scope"] not in ["project","sample"]:
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
            
            # # Not defined yet
            # raise AssertionExcept("Assembly for scope 'project' is not defined yet in %s\n")
            # # See clc_assembl for definition, also in step_sample_initiation() of clc_assembl...
            
                            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Megahit requires that the sample dir not exist! Removing:
            self.script += "rm -rf %s\n\n" % use_dir

            out_prefix = self.sample_data["Title"] + self.file_tag
            
            self.script += self.get_script_const()
 
            f_reads_csl = ""
            r_reads_csl = ""
            s_reads_csl = ""

            f_reads_csl = ",".join([self.sample_data[sample]["fastq.F"] for sample in self.sample_data["samples"] if "fastq.F" in self.sample_data[sample]])
            r_reads_csl = ",".join([self.sample_data[sample]["fastq.R"] for sample in self.sample_data["samples"] if "fastq.R" in self.sample_data[sample]])
            s_reads_csl = ",".join([self.sample_data[sample]["fastq.S"] for sample in self.sample_data["samples"] if "fastq.S" in self.sample_data[sample]])
            # for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            #     if "PE" in self.sample_data[sample]["type"]:
            #         f_reads_csl += "%s,\\\n\t\t" % self.sample_data[sample]["fastq.F"]
            #         r_reads_csl += "%s,\\\n\t\t" % self.sample_data[sample]["fastq.R"]
            #     if "SE" in self.sample_data[sample]["type"]:
            #         s_reads_csl += "%s,\\\n\t\t" % self.sample_data[sample]["fastq.S"]
            #     if "PE" not in self.sample_data[sample]["type"] and "SE" not in self.sample_data[sample]["type"]:
            #         raise AssertionExcept("Strange type configuration for sample\n" ,sample)

            # Interlaced reads to treated here. Maybe one day...
            if f_reads_csl:
                self.script += "-1 " + f_reads_csl.strip(",\\\n\t\t") + " \\\n\t"
                self.script += "-2 " + r_reads_csl.strip(",\\\n\t\t") + " \\\n\t"
            if s_reads_csl:
                self.script += "-r " + s_reads_csl.strip(",\\\n\t\t") + " \\\n\t"

            self.script += "--out-dir %s \\\n\t" % use_dir
            self.script += "--out-prefix %s \n\n" % out_prefix
            
            # Store results to fasta and assembly slots:
            self.sample_data["project_data"]["fasta.nucl"] = self.base_dir + out_prefix + ".contigs.fa"
            self.sample_data["project_data"][self.get_step_step() + "_contigs"] = self.sample_data["project_data"]["fasta.nucl"]
            self.sample_data["project_data"]["assembler"] = self.get_step_step()

            self.stamp_file(self.sample_data["project_data"][self.get_step_step() + "_contigs"])
                

                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                    
        
        else: # self.params["scope"] == "sample"
        
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

                # Megahit requires that the sample dir not exist! Removing:
                self.script += "rm -rf %s\n\n" % use_dir

                out_prefix = sample + self.file_tag
                
                self.script += self.get_script_const()
                self.script += "--out-dir %s \\\n\t" % sample_dir
                self.script += "--out-prefix %s \\\n\t" % out_prefix
     

                if "PE" in self.sample_data[sample]["type"]:
                    self.script += "-1 %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                    self.script += "-2 %s \n\n" % self.sample_data[sample]["fastq.R"]
                elif "SE" in self.sample_data[sample]["type"]:
                    self.script += "-r %s \n\n" % self.sample_data[sample]["fastq.S"]
                elif "PE" in self.sample_data[sample]["type"] and "SE" in self.sample_data[sample]["type"]:       # Mixed!!
                    self.script += "-1 %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                    self.script += "-2 %s \\\n\t" % self.sample_data[sample]["fastq.R"]
                    self.script += "-r %s \n\n" % self.sample_data[sample]["fastq.S"]
                else:
                    raise AssertionExcept("Strange type configuration for sample\n" ,sample)

                # Store results to fasta and assembly slots:
                self.sample_data[sample]["fasta.nucl"] = sample_dir + out_prefix + ".contigs.fa"
                self.sample_data[sample][self.get_step_step() + "_contigs"] = self.sample_data[sample]["fasta.nucl"]
                self.sample_data[sample]["assembler"] = self.get_step_step()

                self.stamp_file(self.sample_data[sample][self.get_step_step() + "_contigs"])
                    
    
                    
                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                
                self.create_low_level_script()
                        

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        pass