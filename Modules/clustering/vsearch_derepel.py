# -*- coding: UTF-8 -*-
""" 
Module ``vsearch_derepel``
--------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running vsearch read dereplication:

Performs dereplication on fastq and fasta files.

..Note:: Dereplication with vsearch is not defined on paired end reads.

    At the moment, this module is defined only for ``fasta.nucl`` or for ``fastq.S``.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files in the following slots:

    * ``sample_data[<sample>]["fastq.S"]``
    
* or fasta files the following slot:

    * ``sample_data[<sample>]["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output **fasta** file in the following slots:

    * ``self.sample_data[<sample>]["fasta.nucl"]``
    * ``self.sample_data[<sample>]["vsearch_derepl"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "sample | project", "Which file to use for dereplication: sample-wise or project-wise files"
    "uc", "-", "Save UCLUST-like dereplication output? (see --uc in manual) "
    "type", "derep_fulllength | derep_prefix", "Type of derelpication strategy. See manual"
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For external index:

::

    derepel_proj:
        module: vsearch_derepel
        base: merge_proj
        script_path: '{Vars.vsearch_path}/vsearch'
        scope: project
        type: derep_fulllength
        uc: 
        redirects:
            --sizein:
            --sizeout:


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rognes, T., Flouri, T., Nichols, B., Quince, C. and Mahé, F., 2016. **VSEARCH: a versatile open source tool for metagenomics**. *PeerJ*, 4, p.e2584.

"""


import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_vsearch_derepel(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
    
        if "scope" not in self.params:
            raise AssertionExcept("You must specify 'scope'\n")
            
        if "type" not in self.params:
            raise AssertionExcept("You must specify 'type': Either 'derep_fulllength' or 'derep_prefix'")
        if self.params["type"] not in ["derep_fulllength" , "derep_prefix"]:
            raise AssertionExcept("You must specify 'type': Either 'derep_fulllength' or 'derep_prefix'")
        
        for type in ["derep_fulllength" , "derep_prefix"]:
            if type in self.params["redir_params"].keys():
                self.write_warning()
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        
        if self.params["scope"] == "sample":
        # Initializing a "fasta" dict for each sample:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Use the fasta->nucl by preference (should be a user defined param...)
                if "fasta.nucl" in self.sample_data[sample] and "fastq.S" in self.sample_data[sample]:
                    if "source" not in self.params:
                        self.write_warning("Both 'nucl' and 'fastq.S' exist and no 'source' param passed. Using 'nucl'")
                        self.params["source"] = "fasta.nucl"
                if "fasta.nucl" not in self.sample_data[sample] and "fastq.S" not in self.sample_data[sample]:
                    raise AssertionExcept("Neither 'nucl' nor 'fastq.S' exist!")
                if "fasta.nucl" in self.sample_data[sample]:
                    self.params["source"] = "fasta.nucl"
                if "fastq.S" in self.sample_data[sample]:
                    self.params["source"] = "fastq.S"

        elif self.params["scope"] == "project":
            # Use the fasta->nucl by preference (should be a user defined param...)
            if "fasta.nucl" in self.sample_data and "fastq.S" in self.sample_data:
                self.write_warning("Both 'fasta.nucl' and 'fastq.S' exist. Using 'nucl'")
                self.params["source"] = "fasta.nucl"
            if "fasta.nucl" not in self.sample_data and "fastq.S" not in self.sample_data:
                raise AssertionExcept("Neither 'fasta.nucl' nor 'fastq.S' exist!")
            if "fasta.nucl" in self.sample_data:
                self.params["source"] = "fasta.nucl"
            if "fastq.S" in self.sample_data:
                self.params["source"] = "fastq.S"


        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        pass
        
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
        if self.params["scope"] == "project":

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            # output_prefix = sample + "_bowtie2_map"

             
            input_file = self.sample_data[self.params["source"]]
            
            output_prefix = os.path.basename(input_file)
            
            # Get constant part of script:
            self.script += self.get_script_const()
            if self.params["type"] == "derep_fulllength":
                self.script += "--derep_fulllength %s \\\n\t" % input_file
            elif self.params["type"] == "derep_prefix":
                self.script += "--derep_prefix %s \\\n\t" % input_file
            if "uc" in self.params:
                self.script += "--uc %s.uc \\\n\t" % (use_dir + output_prefix)
            self.script += "--output %s.fasta \\\n\t" % (use_dir + output_prefix)
            
            

            self.sample_data["fasta.nucl"] = "%s%s.fasta" % (self.base_dir , output_prefix)
            self.sample_data["vsearch_derepl"] = self.sample_data["fasta.nucl"]
            self.stamp_file(self.sample_data["fasta"]["fasta.nucl"])
                    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                            
        else:
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
                # output_prefix = sample + "_bowtie2_map"

                input_file = self.sample_data[sample][self.params["source"]]

                output_prefix = os.path.basename(input_file)

                # Get constant part of script:
                self.script += self.get_script_const()
                if self.params["type"] == "derep_fulllength":
                    self.script += "--derep_fulllength %s \\\n\t" % input_file
                elif self.params["type"] == "derep_prefix":
                    self.script += "--derep_prefix %s \\\n\t" % input_file
                if "uc" in self.params:
                    self.script += "--uc %s.uc \\\n\t" % (use_dir + output_prefix)
                self.script += "--output %s.fasta \\\n\t" % (use_dir + output_prefix)
                
                

                self.sample_data[sample]["fasta.nucl"] = "%s%s.fasta" % (sample_dir , output_prefix)
                self.sample_data[sample]["vsearch_derepl"] = self.sample_data[sample]["fasta.nucl"]
                self.stamp_file(self.sample_data[sample]["fasta.nucl"])
                        
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
           
                
                
                self.create_low_level_script()
                    
