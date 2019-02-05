# -*- coding: UTF-8 -*-
""" 
``fastqc_html`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running fastqc. 

Creates scripts that run fastqc on all available fastq files.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* puts fastqc output files in the following slots:
        
    * ``sample_data[<sample>]["fastqc_fastq.F_html"]``
    * ``sample_data[<sample>]["fastqc_fastq.R_html"]``
    * ``sample_data[<sample>]["fastqc_fastq.S_html"]``
            
* puts fastqc zip files in the following slots:
        
    * ``sample_data[<sample>]["fastqc_fastq.F_zip"]``
    * ``sample_data[<sample>]["fastqc_fastq.R_zip"]``
    * ``sample_data[<sample>]["fastqc_fastq.S_zip"]``
            
 

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::

    fqc_merge1:
        module: fastqc_html
        base: merge1
        script_path: /path/to/FastQC/fastqc
        qsub_params:
            -pe: shared 15
        redirects:
            --threads: 15

 

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Andrews, S., 2010. FastQC: a quality control tool for high throughput sequence data. 

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_fastqc_html(Step):

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "fasqc"

        if "scope" not in self.params:
            self.params["scope"] = "sample"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:      # Getting list of samples out of samples_hash

        # Assert that all samples have reads files:
        # for sample in self.sample_data["samples"]:
            if not [x for x in list(self.sample_data[sample].keys()) if x in ["fastq.F", "fastq.R", "fastq.S"]]:
                raise AssertionExcept("No read files defined\n", sample)

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        if "sum_script" in list(self.params.keys()):
            raise AssertionExcept("Option 'sum_script' no longer supported. Use MultiQC instead!")
            # self.script = "%(script)s \\\n\t-d %(indir)s \\\n\t-o %(outdir)s\n\n" % \
            #         {"script" : self.params["sum_script"],
            #          "indir"  : self.base_dir,
            #          "outdir" : self.base_dir}

        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
        # for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
        #
        #     # Name of specific script:
        #     self.spec_script_name = self.set_spec_script_name(sample)
        #     self.script = ""
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            
            self.script += self.get_script_const()

            self.script += "--outdir " + use_dir
            
            for direction in ("fastq.F","fastq.R","fastq.S"):
                if direction in list(self.sample_data[sample].keys()):
                    self.script += " \\\n\t" + self.sample_data[sample][direction] 
            self.script += "\n\n";
            
            
            # Create temporary dict to store the output file names:
            temp_dict = {}
            for direction in ("fastq.F","fastq.R","fastq.S"):
                if direction in list(self.sample_data[sample].keys()):
                    # temp_dict[direction] = {}
                    file_basename = os.path.basename(self.sample_data[sample][direction])
                    for type in ["zip","html"]:
                        slot_name = "fastqc_{direction}_{type}".format(direction=direction, type=type)
                        self.sample_data[sample][slot_name] = "{dir}{bn}_fastqc.{type}".format(dir=sample_dir,
                                                                                               bn=file_basename,
                                                                                               type=type)

            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir, sample_dir)       # Sees to copying local files to final destination (and other stuff)

            
                    
            
            self.create_low_level_script()
                    
