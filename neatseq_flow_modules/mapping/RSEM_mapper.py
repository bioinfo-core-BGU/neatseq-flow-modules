# -*- coding: UTF-8 -*-
""" 
``RSEM_prep``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``rsem-calculate-expression``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``  (``scope`` = ``project``)
    * ``sample_data[<sample>]["fasta.nucl"]``   (``scope`` = ``sample``)
    
* If neither exists, please supply ``reference`` parameter.
       

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Puts output index files in one of the following slot:

    * ``self.sample_data[<sample>]["RSEM_index"]``
    * ``self.sample_data["RSEM_index"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Where to take the reference from"
    "reference", "path to reference", "Use this fasta file. See the definition for reference_fasta_file(s) in the ARGUMENTS section of rsem-prepare-reference help"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::


    RSEM_prep_ind:
        module:             RSEM_prep
        base:               merge1
        script_path:        /path/to/RSEM
        reference:              /path/to/fasta
        redir_params:
            --gtf:          /path/to/gtf
            --transcript-to-gene-map: /path/to/map_file
    
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_RSEM_mapper(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        if "scope" not in self.params:
            raise AssertionExcept("Please supply a scope parameter: either 'sample' or 'project'!")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:
                if "RSEM_index" not in self.sample_data[sample]:
                    raise AssertionExcept("No RSEM_index exists for RSEM mapper in sample!", sample)
        else:    #if self.params["scope"] == "project":
            if "RSEM_index" not in self.sample_data:
                raise AssertionExcept("No RSEM_index exists for RSEM mapper!")

            
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        # Each iteration must define the following class variables:
            # self.spec_script_name
            # self.script
        
            
    
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define location and prefix for output files:
            output_prefix = sample + "_RSEM"
            # Getting alignment file, if exists, in this order: bam, sam, cram
            alignment = None
            for type in ["bam","sam","cram"]:
                if type in self.sample_data[sample]:
                    alignment = self.sample_data[sample][type]
                 
            # Get constant part of script:
            self.script += self.get_script_const()

            if alignment:
                self.script += "--alignment %s \\\n\t" % alignment
            elif "fastq.F" in self.sample_data[sample]:
                self.script += "--paired-end \\\n\t\t%s \\\n\t\t%s \\\n\t" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
            elif "fastq.S" in self.sample_data[sample]:
                self.script += "%s \\\n\t" % self.sample_data[sample]["fastq.S"]
            else:
                raise AssertionExcept("No fastq file. Thats really really strange...\n")


            if self.params["scope"] == "sample":
                self.script += "%s \\\n\t" % self.sample_data[sample]["RSEM_index"]
            else:    #if self.params["scope"] == "project":
                self.script += "%s \\\n\t" % self.sample_data["RSEM_index"]

            # Savinf bam files:
            if "--output-genome-bam" in self.params["redir_params"].keys():
                if "--sort-bam-by-read-name" in self.params["redir_params"] or "--sort-bam-by-coordinate" in self.params["redir_params"]:
                    self.sample_data[sample]["genome.unsorted.bam"] = sample_dir + sample + ".genome.bam"
                    self.sample_data[sample]["genome.bam"] = sample_dir + sample + ".genome.sorted.bam"
                else:
                    self.sample_data[sample]["genome.bam"] = sample_dir + sample + ".genome.bam"
            if "--sort-bam-by-read-name" in self.params["redir_params"] or "--sort-bam-by-coordinate" in self.params["redir_params"]:
                self.sample_data[sample]["transcript.unsorted.bam"] = sample_dir + sample + ".transcript.bam"
                self.sample_data[sample]["transcript.bam"] = sample_dir + sample + ".transcript.sorted.bam"
            else:
                self.sample_data[sample]["transcript.bam"] = sample_dir + sample + ".transcript.bam"

            
            self.sample_data[sample]["genes.counts"] = sample_dir + sample + ".isoforms.results"

            self.sample_data[sample]["isoforms.counts"] = sample_dir + sample + ".genes.results"
    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
        
