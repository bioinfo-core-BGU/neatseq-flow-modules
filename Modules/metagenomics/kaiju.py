# -*- coding: UTF-8 -*-
""" 
``kaiju``                            
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

    * ``self.sample_data["krona"]``

    
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
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"

class Step_kaiju(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".kaiju.out"


            
        # Checking this once and then applying to each sample:
        try:
            self.params["redir_params"]
        except KeyError:
            raise AssertionExcept("You must specify -t and -f as redirected params.\n")
            
        if "-t" not in self.params["redir_params"].keys() or "-f" not in self.params["redir_params"].keys():
            raise AssertionExcept("You must specify -t and -f as redirected params.\n")


        if "kaiju2krona" in self.params.keys():
            if self.params["kaiju2krona"] == None:
                self.params["kaiju2krona"] = "%s%s%s" % (os.path.basename(self.params["script_path"]), \
                                                            os.sep, \
                                                            "kaiju2krona")
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
 
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.make_sample_file_index()   # see definition below
        
        try:
            self.params["ktImportText_path"]
        except KeyError:
            self.write_warning("You did not supply a ktImportText_path. Will not create krona reports...\n")
            self.script = ""
        else:
            krona_report_fn = self.base_dir + self.sample_data["Title"] + "_krona_report.html"
            self.script = "# Creating krona html reports\n"
            # Adding env and setenv lines to script
            self.script += self.get_setenv_part()
            # Main part of script:
            self.script += "%s \\\n\t" % self.params["ktImportText_path"]
            self.script += "-o %s \\\n\t" % krona_report_fn
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s,%s \\\n\t" % (self.sample_data[sample]["classification"],sample)
            # Removing extra \\
            self.script = self.script.rstrip("\\\n\t") 

            # Storing and stamping results:
            self.sample_data["krona"] = krona_report_fn
            self.stamp_file(self.sample_data["krona"])
            
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
            output_filename = "".join([use_dir , sample , self.file_tag])

            self.script += self.get_script_const()
            
            # Adding reads
            if "PE" in self.sample_data[sample]["type"]:
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                self.script += "-j %s \\\n\t" % self.sample_data[sample]["fastq.R"]
            elif "SE" in self.sample_data[sample]["type"]:
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["fastq.S"]
            else:
                self.write_warning("kaiju is not defined on mixed PE/SE samples. Using only PE data!\n")
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                self.script += "-j %s \\\n\t" % self.sample_data[sample]["fastq.R"]

            self.script += "-o %s\n\n" % output_filename

            # Storing the output file in $samples_hash
            self.sample_data[sample]["raw_classification"]        = "%s" % (sample_dir + os.path.basename(output_filename))
            self.stamp_file(self.sample_data[sample]["raw_classification"])

            if "kaiju2krona" in self.params.keys():
                
                self.script += "# Creating text report for krona\n"
                self.script += "%s \\\n\t" % self.params["kaiju2krona"]
                self.script += "-t %s \\\n\t" % self.params["redir_params"]["-t"]
                self.script += "-n %s \\\n\t" % self.params["names_dmp"]
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["raw_classification"]
                self.script += "-o %s.4krona.txt \n\n" % self.sample_data[sample]["raw_classification"]
                

                self.sample_data[sample]["classification"] = "%s.4krona.txt" % self.sample_data[sample]["raw_classification"]
                self.stamp_file(self.sample_data[sample]["classification"])

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "kaiju_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tkaiju_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["raw_classification"]))
                
        self.sample_data["kaiju_file_index"] = self.base_dir + "kaiju_files_index.txt"
        