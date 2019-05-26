# -*- coding: UTF-8 -*-
""" 
``metaphlan2``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``metaphlan2``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either paired end or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``metaphlan2`` output files in:  

    * ``self.sample_data[<sample>]["raw_classification"]``

* If 
    
* If ``ktImportText_path`` parameter was passed, puts the krona reports in 

    * ``self.sample_data["project_data"]["krona"]``

* If ``merge_metaphlan_tables`` was passed, puts the merged reports in 

    * ``self.sample_data["project_data"]["merged_metaphlan2"]``


* If '--biom' is set in ``redirects``, the biom table is put in:

    * ``self.sample_data[<sample>]["biom_table"]``
    
* If '--bowtie2out' is set in ``redirects``, the SAM file is put in:

    * ``self.sample_data[<sample>]["sam"]``
    
* If  'metaphlan2krona_path' is set:

    * ``self.sample_data[<sample>]["classification"]``

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ktImportText_path",      "", "Path to ktImportText."
    "merge_metaphlan_tables", "", "Path to merge_metaphlan_tables.py. If not specified, will derive it from the location of ``metaphlan2``"
    "metaphlan2krona_path",   "", "Path to metaphlan2krona.py"


    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    metph1:
        module: metaphlan2
        base: trim1
        script_path: {Vars.paths.metaphlan2}
        ktImportText_path: /path/to/ktImportText
        merge_metaphlan_tables: 
        metaphlan2krona_path:   /path/to/metaphlan2krona.py
        redirects:
            --biom: 
            --bowtie2_exe: /path/to/bowtie2
            --bowtie2db: /path/to/database
            --bowtie2out:
            --input_type: fastq
            --mdelim: ';'
            --mpa_pkl: /path/to/mpa_v20_m200.pkl
            
            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Truong, D.T., Franzosa, E.A., Tickle, T.L., Scholz, M., Weingart, G., Pasolli, E., Tett, A., Huttenhower, C. and Segata, N., 2015. **MetaPhlAn2 for enhanced metagenomic taxonomic profiling**. *Nature methods*, 12(10), pp.902-903.

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_metaphlan2(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".metaphlan.out"

        # self.auto_redirs = "--input_type".split(" ")

        if "--input_type" in self.params["redir_params"]:
            self.write_warning("At the moment metaphlan supports only --input_type fastq. Ignoring the value you passed\n")
        
        self.params["redir_params"]["--input_type"] = "fastq"
        
        if "--bowtie2out" in self.params["redir_params"] and self.params["redir_params"]["--bowtie2out"]:
            self.write_warning("Ignoring the value you passed for --bowtie2out.\nWill store data in sample specific location\n")

        if "--biom" in self.params["redir_params"] and self.params["redir_params"]["--biom"]:
            self.write_warning("Ignoring the value you passed for --biom.\nWill store data in sample specific location\n")

        if "merge_metaphlan_tables" in self.params:
            if not isinstance(self.params["merge_metaphlan_tables"], dict):
                if self.params["merge_metaphlan_tables"]:
                    raise AssertionExcept("'merge_metaphlan_tables' must be empty or a block with 'path' and optionally 'redirects'")
                else:
                    self.params["merge_metaphlan_tables"] = dict()
            if "path" not in self.params["merge_metaphlan_tables"]:
                self.params["merge_metaphlan_tables"]["path"] = os.sep.join([os.path.basename(self.params["script_path"]),
                                                                             "utils/merge_metaphlan_tables.py"])
                self.write_warning("You did not provided a path in 'merge_metaphlan_tables'. Using '{path}'".
                                   format(path=self.params["merge_metaphlan_tables"]["path"]))

        if "ktImportText" in self.params:
            if not isinstance(self.params["ktImportText"], dict) or "path" not in self.params["ktImportText"]:
                raise AssertionExcept("Please include a 'path' in the 'ktImportText' block.")
        else:
            self.write_warning("You did not supply a 'ktImportText' block. Will not create krona reports...\n")


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
 
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        

        try:
            self.params["ktImportText"]
        except KeyError:
            self.write_warning("You did not supply a 'ktImportText' block. Will not create krona reports...\n")
            self.script = ""
        else:
            if "redirects" in self.params["ktImportText"]:
                redirects = " \\\n\t" + " \\\n\t".join(
                    [key + " " + (val if val else "")
                     for key, val
                     in self.params["ktImportText"]["redirects"].items()])
            else:
                redirects = ""

            self.make_sample_file_index()  # see definition below

            krona_report_fn = self.base_dir + self.sample_data["Title"] + "_krona_report.html"

            # Adding env to script
            self.script = self.get_setenv_part()
            self.script += "# Creating krona html reports\n"
            # if "env" in self.params:
            #     self.script += "env %s \\\n\t" % self.params["env"]
            # Main part of script:
            self.script += """
# Executing ktImportText on metaphlan2 results:
{path} \\{redirs}
\t-o {output}
""".format(path=self.params["ktImportText"]["path"],
           redirs=redirects,
           output=krona_report_fn)

            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s,%s \\\n\t" % (self.sample_data[sample]["classification"], sample)
            self.script = self.script.rstrip("\\\n\t") 
            self.script += "\n\n\n"
            # Storing and stamping results:
            self.sample_data["project_data"]["krona"] = krona_report_fn
            self.stamp_file(self.sample_data["project_data"]["krona"])
            
        if "merge_metaphlan_tables" in self.params:

            if "redirects" in self.params["merge_metaphlan_tables"]:
                redirects = " \\\n\t" + " \\\n\t".join(
                    [key + " " + (val if val else "")
                     for key, val
                     in self.params["merge_metaphlan_tables"]["redirects"].items()])
            else:
                redirects = ""

            self.script += """# Merging all metaphlan2 reports into single table
{path} \\{redirs}
\t{files} \\
\t> {dir}{outfn}""".format(path=self.params["merge_metaphlan_tables"]["path"],
                           redirs=redirects,
                           files="\\\n\t".join([self.sample_data[sample]["raw_classification"]
                                                for sample
                                                in self.sample_data["samples"]]),
                           outfn=self.sample_data["Title"] + "_merged_table.txt",
                           dir= self.base_dir)

            # for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            #     self.script += "%s \\\n\t" % self.sample_data[sample]["raw_classification"]
            # self.script += "> %s%s \n\n" %

            self.sample_data["project_data"]["merged_metaphlan2"] = (self.base_dir, self.sample_data["Title"]+"_merged_table.txt")
        
        
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
    
        
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
                
            # Define output filename
            output_filename = "".join([use_dir , sample , self.file_tag])

            # If user passed bowtie2out or biom, change the value to sample specific values:
            # These are then added with redirected params
            if "--bowtie2out" in list(self.params["redir_params"].keys()):
                self.params["redir_params"]["--bowtie2out"] = "%s.sam" % output_filename

            if "--biom" in list(self.params["redir_params"].keys()):
                self.params["redir_params"]["--biom"] = "%s.biom" % output_filename

            self.script += self.get_script_const()
            self.script += "--sample_id %s \\\n\t" % sample
            
            # Adding reads

            if "fastq.F" in self.sample_data[sample]:
                self.script += "%s,%s\\\n\t" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
            elif "fastq.S" in self.sample_data[sample]:
                self.script += "%s \n\n" % self.sample_data[sample]["fastq.S"]
            else:
                raise AssertionExcept("Weird cmobination of reads!\n")

            self.script += "%s\n\n" % output_filename

            # Storing the output file in $samples_hash
            self.sample_data[sample]["raw_classification"]        = "%s" % (sample_dir + os.path.basename(output_filename))
            self.stamp_file(self.sample_data[sample]["raw_classification"])

            # Storing and stamping biom and sam files if requested
            if "--biom" in list(self.params["redir_params"].keys()):
                # if "qiime" not in self.sample_data[sample].keys():
                    # self.sample_data[sample]["qiime"] = dict()
                self.sample_data[sample]["biom_table"] = "%s.biom" % (sample_dir + os.path.basename(output_filename))
                self.stamp_file(self.sample_data[sample]["biom_table"])
            if "--bowtie2out" in list(self.params["redir_params"].keys()):
                # if "mapping" not in self.sample_data[sample].keys():
                    # self.sample_data[sample]["mapping"] = dict()
                self.sample_data[sample]["sam"] = "%s.sam" % (sample_dir + os.path.basename(output_filename))
                self.stamp_file(self.sample_data[sample]["sam"])

            if "metaphlan2krona" in self.params:
                
                self.script += "# Creating text report for krona\n"
                self.script += "%s \\\n\t" % self.params["metaphlan2krona"]["path"]
                self.script += "-p %s \\\n\t" % self.sample_data[sample]["raw_classification"]
                self.script += "-k %s.4krona.txt\n\n" % self.sample_data[sample]["raw_classification"]

                self.sample_data[sample]["classification"] = "%s.4krona.txt" % (sample_dir + os.path.basename(output_filename))
                self.stamp_file(self.sample_data[sample]["classification"])

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "metaphlan2_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tmetaphlan2_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["classification"]))
                
        self.sample_data["project_data"]["metaphlan2_files_index"] = self.base_dir + "metaphlan2_files_index.txt"