# -*- coding: UTF-8 -*-
""" 
``mash_dist`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

    
Requires:
~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``

* OR fastq files in one of the following slots (merge fastq files first with mash_sketch or otherwise):

    * ``sample_data[<sample>]["fastq"]``
    * ``sample_data["fastq"]``

* OR sketch files in one of the following slots:

    * ``sample_data[<sample>]["msh.fastq"]``
    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data["msh.fastq"]``
    * ``sample_data["msh.fasta"]``
    
    
Output:
~~~~~~~~~~~~~

* puts 'msh' output files in the following slots for (scope=sample):

    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data[<sample>]["msh.fastq"]``

* puts 'msh' output files in the following slots for (scope=project and scope=all_samples):

    * ``sample_data[<sample>]["mash.dist.table"]``
    * ``sample_data["mash.dist.table"]``
    



Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10
    
    "reference", "", "A block including 'path' or 'scope', 'type' and optionally 'msh'"
    "query", "", "A block including 'scope' (sample, project or all_samples), 'type' and optionally 'msh'"


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::
    dist: 
        module:         mash_dist
        base:           [sketch_proj,sketch_smp]
        script_path:    "{Vars.paths.mash} dist"
        reference:
            # path:   /path/to/ref1
            type:           fastq
            msh:
        query:
            scope:          all_samples
            type:           fastq
            msh:


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step, AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_mash_dist(Step):
    
    auto_redirs = "-l".split(" ")

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".mash.tbl"


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        ## TODO: Check that files exist for the analysis
        pass
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    
    def build_scripts(self):
        """ This is the actual script building function
            
        """


          
        if self.params["query"]["scope"] in ["all_samples","project"]:
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()
        

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
        """
   
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
            
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if "path" in self.params["reference"]:
                ref_path = self.params["reference"]["path"]
            else:
                if "msh" in self.params["reference"]:
                    if "type" in self.params["reference"] and self.params["reference"]["type"]=="fasta":
                        try:
                            ref_path = self.sample_data["msh.fasta"]
                        except KeyError:
                            raise AssertionExcept("No project fasta mash sketch")
                    else:
                        try:
                            ref_path = self.sample_data["msh.fastq"]
                        except KeyError:
                            raise AssertionExcept("No project fastq mash sketch")

                else:
                    if "type" in self.params["reference"] and self.params["reference"]["type"]=="fasta":
                        try:
                            ref_path = self.sample_data["fasta.nucl"]
                        except KeyError:
                            raise AssertionExcept("No project fasta.nucl file")
                    else:
                        try:
                            ref_path = self.sample_data["fastq"]
                        except KeyError:
                            raise AssertionExcept("No project 'fastq' file")
            ################ Setting query path
            if "msh" in self.params["query"]:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = self.sample_data[sample]["msh.fasta"]
                    except KeyError:
                        raise AssertionExcept("No sample fasta mash sketch",sample)
                else:
                    try:
                        query_path = self.sample_data[sample]["msh.fastq"]
                    except KeyError:
                        raise AssertionExcept("No sample fastq mash sketch",sample)

            else:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = self.sample_data[sample]["fasta.nucl"]
                    except KeyError:
                        raise AssertionExcept("No sample fasta.nucl file",sample)
                else:
                    try:
                        query_path = self.sample_data[sample]["fastq"]
                    except KeyError:
                        raise AssertionExcept("No sample 'fastq' file",sample)

            print query_path
            print ref_path

            
            output_filename = sample + self.file_tag
            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            self.script += self.get_script_const()
            self.script += "%s \\\n\t" % ref_path
            self.script += "%s \\\n\t" % query_path
            self.script += "> %s \n\n" % (use_dir + output_filename)
            
                        
            # Store results table
            self.sample_data[sample]["mash.dist.table"] = (sample_dir + output_filename)
            self.stamp_file(self.sample_data[sample]["mash.dist.table"])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
            
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
        self.script = ""

        if "path" in self.params["reference"]:
            ref_path = self.params["reference"]["path"]
        else:
            if "msh" in self.params["reference"]:
                if "type" in self.params["reference"] and self.params["reference"]["type"]=="fasta":
                    try:
                        ref_path = self.sample_data["msh.fasta"]
                    except KeyError:
                        raise AssertionExcept("No project fasta mash sketch")
                else:
                    try:
                        ref_path = self.sample_data["msh.fastq"]
                    except KeyError:
                        raise AssertionExcept("No project fastq mash sketch")

            else:
                if "type" in self.params["reference"] and self.params["reference"]["type"]=="fasta":
                    try:
                        ref_path = self.sample_data["fasta.nucl"]
                    except KeyError:
                        raise AssertionExcept("No project fasta.nucl file")
                else:
                    try:
                        ref_path = self.sample_data["fastq"]
                    except KeyError:
                        raise AssertionExcept("No project 'fastq' file")
        ################ Setting query path
        if self.params["query"]["scope"] == "project":
            if "msh" in self.params["query"]:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = self.sample_data["msh.fasta"]
                    except KeyError:
                        raise AssertionExcept("No project fasta mash sketch")
                else:
                    try:
                        query_path = self.sample_data["msh.fastq"]
                    except KeyError:
                        raise AssertionExcept("No project fastq mash sketch")

            else:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = self.sample_data["fasta.nucl"]
                    except KeyError:
                        raise AssertionExcept("No project fasta.nucl file")
                else:
                    try:
                        query_path = self.sample_data["fastq"]
                    except KeyError:
                        raise AssertionExcept("No project 'fastq' file")
                        
        else:  # scope = "all_samples"
            if "msh" in self.params["query"]:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = " ".join([self.sample_data[sample]["msh.fasta"] for sample in self.sample_data["samples"]])
                    except KeyError:
                        raise AssertionExcept("A sample is missing a fasta mash sketch")
                else:
                    try:
                        query_path = " ".join([self.sample_data[sample]["msh.fastq"] for sample in self.sample_data["samples"]])

                    except KeyError:
                        raise AssertionExcept("A sample is missing a fastq mash sketch")

            else:
                if "type" in self.params["query"] and self.params["query"]["type"]=="fasta":
                    try:
                        query_path = " ".join([self.sample_data[sample]["fasta.nucl"] for sample in self.sample_data["samples"]])
                    except KeyError:
                        raise AssertionExcept("A sample is missing a fasta.nucl file")
                else:
                    try:
                        query_path = " ".join([self.sample_data[sample]["fastq"] for sample in self.sample_data["samples"]])

                    except KeyError:
                        raise AssertionExcept("A sample is missing a 'fastq' file")

        output_filename = self.sample_data["Title"] + self.file_tag

        self.script += self.get_script_const()
        self.script += "%s \\\n\t" % ref_path
        self.script += "%s \\\n\t" % query_path
        self.script += "> %s \n\n" % (use_dir + output_filename)
        
                    
        # Store results table
        self.sample_data["mash.dist.table"] = (self.base_dir + output_filename)
        self.stamp_file(self.sample_data["mash.dist.table"])


        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        self.create_low_level_script()
    
