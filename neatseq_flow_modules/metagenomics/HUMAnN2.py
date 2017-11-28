# -*- coding: UTF-8 -*-
""" 
``HUMAnN2``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``HUMAnN2``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either forward or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


self.sample_data[sample]["HUMAnN2.genefamilies.norm"]
self.sample_data[sample]["HUMAnN2.pathabundance.norm"]


* Puts the ``HUMAnN2`` output files in:  

    * ``self.sample_data[sample]["HUMAnN2.genefamilies"]``
    * ``self.sample_data[sample]["HUMAnN2.pathabundance"]``
    * ``self.sample_data[sample]["HUMAnN2.pathcoverage"]``

* If "renorm_table" is set in params:

    * ``self.sample_data[sample]["HUMAnN2.genefamilies.norm"]``
    * ``self.sample_data[sample]["HUMAnN2.pathabundance.norm"]``
    
* If "join_tables" is set in params:
    
    * ``self.sample_data["HUMAnN2.genefamilies"]``
    * ``self.sample_data["HUMAnN2.pathabundance"]``
    * ``self.sample_data["HUMAnN2.pathcoverage"]``

* If "join_tables" and "renorm_table" are set in params:
    
    * ``self.sample_data["HUMAnN2.genefamilies.norm"]``
    * ``self.sample_data["HUMAnN2.pathabundance.norm"]``
    * ``self.sample_data["HUMAnN2.pathcoverage"]``
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "renorm_table",      "Empty or parameters to pass to ``humann2_renorm_table``", "Runs ``humann2_renorm_table`` on HUMAnN2 outputs with the specified parameters"
    "join_tables", "", "Runs ``humann2_join_tables`` to gather all sample tables."


    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    HUMAnN2_1:
        module: HUMAnN2
        base: trim1
        script_path: {Vars.paths.humann2}
        join_tables: 
        renorm_table: --units cpm -p
        redirects:
            --bowtie2: /path/to/bowtie2
            --gap-fill: on
            --input-format: fastq
            --metaphlan: {Vars.paths.metaphlan2}
            --minpath: on
            --nucleotide-database: {Vars.databases.chocophlan}
            --protein-database: {Vars.databases.uniref}
            --threads: 30
            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`HUMAnN2 home page <http://huttenhower.sph.harvard.edu/humann2>`_


"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_HUMAnN2(Step):

    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "HUMAnN2"

            
        # # Checking this once and then applying to each sample:
        # try:
            # self.params["redir_params"]
        # except KeyError:
            # raise AssertionExcept("You must specify --mpa_pkl and --bowtie2db as redirected params.\n")
            
        # if "--mpa_pkl" not in self.params["redir_params"].keys() or "--bowtie2db" not in self.params["redir_params"].keys():
            # raise AssertionExcept("You must specify --mpa_pkl and --bowtie2db as redirected params.\n")

        if "--input-format" in self.params["redir_params"].keys():
            self.write_warning("At the moment metaphlan supports only --input-format fastq. Ignoring the value you passed\n")
        
        self.params["redir_params"]["--input-format"] = "fastq"
        
        if "--o-log" in self.params["redir_params"].keys():
            self.write_warning("Ignoring the value you passed for --o-log.\nWill store data in sample specific location\n")

           

        # if "merge_metaphlan_tables" in self.params.keys():
            # if self.params["merge_metaphlan_tables"] == None:
                # self.params["merge_metaphlan_tables"] = "%s%s%s" % (os.path.basename(self.params["script_path"]), \
                                                            # os.sep, \
                                                            # "utils/merge_metaphlan_tables.py")
        
            
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
 
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        if "join_tables" in self.params.keys():
            
            # Get location of humann2 scripts:
            humann2_dir,main_script = os.path.split(self.params["script_path"])

            # Merging genefamilies:
            self.script += "%s \\\n\t" % (humann2_dir + os.sep + "humann2_join_tables")
            self.script += "-i %s \\\n\t" % self.base_dir
            self.script += "--search-subdirectories \\\n\t" 
            if "renorm_table" in self.params.keys():
                output_filename = "genefamilies.norm"
            else:
                output_filename = "genefamilies"
            self.script += "--file_name %s \\\n\t" % output_filename

            self.script += "-o %smerged.%s.tsv \n\n" % (self.base_dir, output_filename)

            ## Storing in dict and stamping
            self.sample_data["HUMAnN2." + output_filename] = "%smerged.%s.tsv" % (self.base_dir,  output_filename)
            self.stamp_file(self.sample_data["HUMAnN2." + output_filename])

            # Merging pathabundance:
            self.script += "%s \\\n\t" % (humann2_dir + os.sep + "humann2_join_tables")
            self.script += "-i %s \\\n\t" % self.base_dir
            self.script += "--search-subdirectories \\\n\t" 
            if "renorm_table" in self.params.keys():
                output_filename = "pathabundance.norm"
            else:
                output_filename = "pathabundance"
            self.script += "--file_name %s \\\n\t" % output_filename

            self.script += "-o %smerged.%s.tsv \n\n" % (self.base_dir, output_filename)

            ## Storing in dict and stamping
            self.sample_data["HUMAnN2." + output_filename] = "%smerged.%s.tsv" % (self.base_dir,  output_filename)
            self.stamp_file(self.sample_data["HUMAnN2." + output_filename])

            # Merging pathcoverage:
            self.script += "%s \\\n\t" % (humann2_dir + os.sep + "humann2_join_tables")
            self.script += "-i %s \\\n\t" % self.base_dir
            self.script += "--search-subdirectories \\\n\t" 
            output_filename = "pathcoverage"
            self.script += "--file_name %s \\\n\t" % output_filename

            self.script += "-o %smerged.%s.tsv \n\n" % (self.base_dir, output_filename)

            ## Storing in dict and stamping
            self.sample_data["HUMAnN2." + output_filename] = "%smerged.%s.tsv" % (self.base_dir,  output_filename)
            self.stamp_file(self.sample_data["HUMAnN2." + output_filename])

            
        
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
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
 
         
            # Define output filename 
            output_filename = "%s_%s" % (sample , self.file_tag)

            # If user passed bowtie2out or biom, change the value to sample specific values:
            # These are then added with redirected params
            if "--o-log" in self.params["redir_params"].keys():
                self.params["redir_params"]["--o-log"] = "%s.log" % output_filename


            self.script += self.get_script_const()
            
            # Adding reads
            if "PE" in self.sample_data[sample]["type"]:
                self.write_warning("PE not defined in HUMAnN2. Using only forward reads.\n")
                self.script += "--input %s \\\n\t" % (self.sample_data[sample]["fastq.F"])
            elif "SE" in self.sample_data[sample]["type"]:
                self.script += "--input %s \n\n" % self.sample_data[sample]["fastq.S"]
            else:
                self.write_warning("metaphlan2 on mixed PE/SE samples is not defined. Using only forward data!\n")
                self.script += "--input %s \\\n\t" % self.sample_data[sample]["fastq.F"]

            self.script += "--output %s \\\n\t" % use_dir
            self.script += "--output-basename %s\n\n" % output_filename

            
            self.sample_data[sample]["HUMAnN2.genefamilies"] = "%s_genefamilies.tsv" % (sample_dir + output_filename)
            self.sample_data[sample]["HUMAnN2.pathabundance"] = "%s_pathabundance.tsv" % (sample_dir + output_filename)
            self.sample_data[sample]["HUMAnN2.pathcoverage"] = "%s_pathcoverage.tsv" % (sample_dir + output_filename)
            
            self.stamp_file(self.sample_data[sample]["HUMAnN2.genefamilies"])
            self.stamp_file(self.sample_data[sample]["HUMAnN2.pathabundance"])
            self.stamp_file(self.sample_data[sample]["HUMAnN2.pathcoverage"])

            
            # Adding code for normalization if required
            
            if "renorm_table" in self.params.keys():
                self.script += "# Adding code for normalizing genefamilies and pathabundance tables\n\n"
                # Get location of humann2 scripts:
                humann2_dir,main_script = os.path.split(self.params["script_path"])

                self.script += "%s \\\n\t" % (humann2_dir + os.sep + "humann2_renorm_table")
                if self.params["renorm_table"]: # If user passed parameters to renorm_table, add them
                    self.script += "%s \\\n\t" % self.params["renorm_table"]
                self.script += "-i %s \\\n\t" % "%s_genefamilies.tsv" % (use_dir + output_filename)
                self.script += "-o %s \n\n" % "%s_genefamilies.norm.tsv" % (use_dir + output_filename)

                self.script += "%s \\\n\t" % (humann2_dir + os.sep + "humann2_renorm_table")
                if self.params["renorm_table"]: # If user passed parameters to renorm_table, add them
                    self.script += "%s \\\n\t" % self.params["renorm_table"]
                self.script += "-i %s \\\n\t" % "%s_pathabundance.tsv" % (use_dir + output_filename)
                self.script += "-o %s \n\n" % "%s_pathabundance.norm.tsv" % (use_dir + output_filename)

                self.sample_data[sample]["HUMAnN2.genefamilies.norm"] = "%s_genefamilies.norm.tsv" % (sample_dir + output_filename)
                self.sample_data[sample]["HUMAnN2.pathabundance.norm"] = "%s_pathabundance.norm.tsv" % (sample_dir + output_filename)
                self.stamp_file(self.sample_data[sample]["HUMAnN2.genefamilies.norm"])
                self.stamp_file(self.sample_data[sample]["HUMAnN2.pathabundance.norm"])
 

            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "HUMAnN2_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tHUMAnN2_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["classification"]))
                
        self.sample_data["HUMAnN2.files_index"] = self.base_dir + "HUMAnN2_files_index.txt"
        