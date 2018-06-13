# -*- coding: UTF-8 -*-
""" 
``mash_sketch`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

    
Requires:
~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``

* or fastq files in the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
* For ``scope = project``, uses project-wide files.
    
Output:
~~~~~~~~~~~~~

* puts 'msh' output files in the following slots for (scope=sample):

    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data[<sample>]["msh.fastq"]``

* puts 'msh' output files in the following slots for (scope=project):

    * ``sample_data["msh.fasta"]``
    * ``sample_data["msh.fastq"]``



Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10
    
    "scope", "project|sample", ""
    "type", "nucl|prot", ""


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_mash_sketch(Step):
    
    auto_redirs = "-o -l".split(" ")

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".msh"

        if self.params["scope"] not in ["sample","project"]:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        if "src_scope" in self.params:
            if self.params["src_scope"] not in ["sample","project"]:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            self.params["src_scope"] = self.params["scope"]

        if "type" not in self.params:
            self.params["type"] = ["fastq","fasta"]
        else:
            if isinstance(self.params["type"], str):
                self.params["type"] = [self.params["type"]]

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        ## TODO: Check that files exist for the analysis
        pass
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.script = ""
        if self.params["scope"]=="project" and self.params["src_scope"] == "sample":
            for type in list(set(self.params["type"])&set(["fastq","fasta"])):

                # Create script only if there are files in files4mashing_lists
                if self.files4mashing_lists[type]:
                    # Define output filename 
                    output_filename = "".join([self.base_dir , self.sample_data["Title"] , ".", type])  # The 'msh' tag get added automatically by mash

                    self.script += self.get_script_const()
                    self.script += "-o %s \\\n\t" % output_filename
                    self.script += "%s \n\n" % " \\\n\t".join(self.files4mashing_lists[type])
                    
                    if "rm_merged" in self.params:  # and len(type_list)>1:
                        self.script += "# Removing temporary merged files\n"
                        for sample in self.files2remove:
                            for file in self.files2remove[sample]:
                                self.script += "rm -rf %s \n\n" % file
                            
                        
                    # Store msh  file:
                    self.sample_data["msh." + type] = (output_filename + ".msh")
                    self.stamp_file(self.sample_data["msh." + type])
                    

    
    def build_scripts(self):
        """ This is the actual script building function
            
        """


          
        if self.params["src_scope"]=="project":
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()
        

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
        """
   
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
            
        self.files4mashing_lists=dict()
        self.files4mashing_lists["fasta"] = list()
        self.files4mashing_lists["fastq"] = list()
        self.files2remove = dict()

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            self.files2remove[sample] = list()

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
            fastq_list = list(set(self.sample_data[sample].keys()) & set(["fastq.F","fastq.R","fastq.S"]))
            fasta_list = list(set(self.sample_data[sample].keys()) & set(["fasta.nucl"]))

            # if fasta_list:  # There are fasta files
                # pass
                
            # if fastq_list:  # There are fastq files

            for type in list(set(self.params["type"])&set(["fastq","fasta"])):
            #list(set(self.params["type"])&set(["fastq","fasta"])):
                if type=="fastq":
                    type_list = fastq_list
                else:
                    type_list = fasta_list
                
                if not type_list:
                    continue
                    
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample,type])
                self.script = ""

                input_filename = ""
                if len(type_list)==1:  # Only one file
                    input_filename = self.sample_data[sample][type_list[0]]
                    final_input_filename = input_filename
                    # self.script += "# Single file in type %s. Not doing concatenation\n\n" % type
                else:
                    input_filename       = use_dir    + sample + ".merged." + type
                    final_input_filename = sample_dir + sample + ".merged." + type
                    self.script += "cat \\\n\t"
                    for spec_type in type_list:
                        self.script += "%s \\\n\t" % self.sample_data[sample][spec_type]
                    self.script += "> %s\n\n" % input_filename
                    self.files2remove[sample].append(input_filename)
                    # This is done in here because sometimes the script is not finalized and 
                    # stamped files get carried from one script to the next.
                    self.stamp_file(final_input_filename)

                self.sample_data[sample][type] = final_input_filename
                

                if self.params["scope"] == "project":
                    self.files4mashing_lists[type].append(final_input_filename)
                    
                if self.params["scope"] == "sample":
                    # Define output filename 
                    output_filename = "%s.%s" % (sample , type)  # The 'msh' tag get added automatically by mash

                    self.script += self.get_script_const()
                    self.script += "-o %s \\\n\t" % (use_dir + output_filename)
                    self.script += "%s \n\n" % input_filename
                    
                    if "rm_merged" in self.params:  # and len(type_list)>1:
                        for file in self.files2remove[sample]:
                            self.script += "rm -rf %s \n\n" % file
                    
                        
                    # Store msh  file:
                    self.sample_data[sample]["msh." + type] = (sample_dir + output_filename + ".msh")
                    self.stamp_file(self.sample_data[sample]["msh." + type])

                if self.script:
                    # Wrapping up function. Leave these lines at the end of every iteration:
                    self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

                    self.create_low_level_script()
            
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        
        type_counter = 0 # Counts the number of types for which info was found. Fails if no files are found

        for type in list(set(self.params["type"])&set(["fastq","fasta"])):
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"],type])
            self.script = ""

            

            if type=="fastq":
                type_list = list(set(self.sample_data.keys()) & set(["fastq.F","fastq.R","fastq.S"]))
            else:
                type_list = list(set(self.sample_data.keys()) & set(["fasta.nucl"]))
            import pdb; pdb.set_trace()

            if not type_list:
                continue
            # Count the number of types for which info was found. Fails if no files are found
            type_counter += 1

            input_filename = ""
            if len(type_list)==1:  # Only one file
                input_filename = self.sample_data[type_list[0]]
                final_input_filename = input_filename
            else:
                input_filename       = use_dir       + self.sample_data["Title"] + ".merged." + type
                final_input_filename = self.base_dir + self.sample_data["Title"] + ".merged." + type
                self.script += "cat \\\n\t"
                for spec_type in type_list:
                    self.script += "%s \\\n\t" % self.sample_data[spec_type]
                self.script += "> %s\n\n" % input_filename
                
                
            # Define output filename 
            output_filename = "%s.%s" % (self.sample_data["Title"] , type)  # The 'msh' tag get added automatically by mash

            self.script += self.get_script_const()
            self.script += "-o %s \\\n\t" % (use_dir + output_filename)
            self.script += "%s \n\n" % input_filename
            
            if "rm_merged" in self.params and len(type_list)>1:
                self.script += "rm -rf %s \n\n" % input_filename
            

            # Store msh  file:
            self.sample_data["msh." + type] = (self.base_dir + output_filename + ".msh")
            self.sample_data[type] = final_input_filename
            self.stamp_file(self.sample_data["msh." + type])
            self.stamp_file(self.sample_data[type])
            

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            self.create_low_level_script()
        

        if type_counter == 0:
            raise AssertionExcept("No source data was found for type %s. Check definitions of 'scope' and 'src_scope'" % ", ".join(list(set([self.params["type"]])&set(["fastq","fasta"]))))
