# -*- coding: UTF-8 -*-
""" 
``kaiju2table``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad


A module for running ``kaiju``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either paired end or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``kaiju`` output files in:  

    * ``self.sample_data[<sample>]["raw_classification"]``
    
* If  'kaiju2krona' is set:

    * ``self.sample_data[<sample>]["classification"]``

* If ``ktImportText_path`` parameter was passed, puts the krona reports in 

    * ``self.sample_data["project_data"]["krona"]``

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ktImportText_path",      "", "Path to ktImportText."
    "kaiju2krona", "", "Path to kaiju2krona. If not specified, will derive it from the ``script_path``"


    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    kaiju1:
        module: kaiju
        base: trim1
        script_path: {Vars.paths.kaiju}
        kaiju2krona: 
        ktImportText_path: {Vars.paths.ktImportText}
        names_dmp: /path/to/kaijudb/names.dmp
        redirects:
            -f: /path/to/kaijudb/kaiju_db.fmi
            -t: /path/to/kaijudb/nodes.dmp
            -z: 40            
            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Menzel, P., Ng, K.L. and Krogh, A., 2016. **Fast and sensitive taxonomic classification for metagenomics with Kaiju**. Nature communications, 7.

"""




import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"



class Step_kaiju2table(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".kaiju.out"

        if "-t" not in self.params["redir_params"]:
            raise AssertionExcept("Please supply Name of nodes.dmp file via '-t' argument (in redirects)")
        if "-n" not in self.params["redir_params"]:
            raise AssertionExcept("Please supply Name of names.dmp file via '-n' argument (in redirects)")

        if "-r" in self.params["redir_params"]:
            if isinstance(self.params["redir_params"]["-r"],list):
                self.levels = self.params["redir_params"]["-r"]
            elif isinstance(self.params["redir_params"]["-r"], str):
                self.levels = re.split("\s*,\s*", self.params["redir_params"]["-r"])
            else:
                raise AssertionExcept("Unknown format of '-r' redirects. Must be either string or list")
            self.params["redir_params"].pop("-r")
        else:
            self.levels = ["phylum", "class", "order", "family", "genus", "species"]


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        # sample_list = self.sample_data["samples"]
        for sample in sample_list:  # Getting list of samples out of samples_hash

            if "raw_classification" not in self.sample_data[sample]:
                raise AssertionExcept("No classification data in sample. Make sure you have a kaiju instance before this instance", sample)


    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """

        if "join_script" in self.params:
            self.script = ""
            print(self.levels)
            for tax_level in self.levels:
                # Define output filename
                print(tax_level)
                output_filename = "{project}.kaiju.summary.{level}.tsv".format(project=self.sample_data["Title"],
                                                                               level=tax_level)

                self.script += """
            
perl {path} \\
\t{samples} \\
\t{files} \\
\t{outfn}

            """.format(path=self.params["join_script"],
                       samples=",".join(self.samples_data["samples"]),
                       files=",".join([self.sample_data[sample]["kaiju.report."+tax_level] for sample in self.samples_data["samples"]]),
                       outfn=self.base_dir+output_filename)
                self.sample_data["project_data"]["kaiju.report."+tax_level] = self.base_dir+output_filename
                # self.stamp_file(self.sample_data["project_data"]["kaiju.report."+tax_level])

    def build_scripts(self):
        """
        """
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        # sample_list = self.sample_data["samples"]
        for sample in sample_list:  # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
 
         
            for tax_level in self.levels:
                # Define output filename
                output_filename = "{sample}.kaiju.summary.{level}.tsv".format(sample=sample, level=tax_level)


                self.script += """
{preamb} \\
\t-r {tax_lev} \\
\t-o {outfn} \\
\t{infn} 


""".format(preamb=self.get_script_const().rstrip(" \\\n\t"),
           tax_lev=tax_level,
           outfn=use_dir+output_filename,
           infn=self.sample_data[sample]["raw_classification"])

                self.sample_data[sample]["kaiju.report."+tax_level] = sample_dir+output_filename
                self.stamp_file(self.sample_data[sample]["kaiju.report."+tax_level])

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        pass