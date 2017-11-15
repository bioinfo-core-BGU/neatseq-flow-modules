# -*- coding: UTF-8 -*-
""" 
``CARD_RGI``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad


A module for running CARD RGI:

RGI is executed on the contigs stored in a Nucleotide fasta file.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A nucleotide fasta file in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* If ``scope`` is set to ``sample``:

    * Puts output files in:
    
        ``sample_data[<sample>]["CARD_RGI.json"]``
        ``sample_data[<sample>]["CARD_RGI.tsv"]``

    * Puts index of output files in:
    
        ``self.sample_data["CARD_RGI.files_index"]``
        
    * If ``merge_script_path`` is specified in parameters, puts the merged file in 
    
        ``self.sample_data["CARD_RGI.merged_reports"]``
        
        
* If ``scope`` is set to ``project``:

    * Puts output files in:
    
        ``sample_data["CARD_RGI.json"]``
        ``sample_data["CARD_RGI.tsv"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "JSON2tsv_script",  "path", "The path to the CARD script for converting the JSON output to tsv 
                                (find 'convertJsonToTSV.py' in your RGI installation)"
    "merge_script_path", "path", "Path to a script that takes an index of RGI output files ('--ind') and a place to put 
                                  the output (--output). This script will be executed in the wrapping up stage.
                                  (Note, the script can take more parameters. These should be passed with the path in 
                                  the parameter files, e.g. 'python /path/to/script --param1 val1 --param2 val2')
                                  If the parameters is not passed, no action will be taken on the output files."
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    rgi_inst:
        module: CARD_RGI
        base: spades1
        script_path: python /path/to/rgi.py
        qsub_params:
            -pe: shared 15
        JSON2tsv_script: python /path/to/convertJsonToTSV.py
        merge_script_path: Rscript /path/to/merge_reports.R --variable bit_score
        orf_to_use: -x
        scope: sample
        redirects:
            -n: 20
            -x: 1

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

McArthur, A.G., Waglechner, N., Nizam, F., Yan, A., Azad, M.A., Baylay, A.J., Bhullar, K., Canova, M.J., De Pascale, G., Ejim, L. and Kalan, L., 2013. **The comprehensive antibiotic resistance database**. *Antimicrobial agents and chemotherapy*, 57(7), pp.3348-3357.

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal_Gordon"

class Step_CARD_RGI(Step):

    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances

        # Making sure only -x or --orf were passed 
        if "-x" in self.params["redir_params"] and "--orf" in self.params["redir_params"]:
            raise AssertionExcept("Please set either -x or --orf (not both)\n")
        
        if "-x" in self.params["redir_params"]:
            self.params["orf_to_use"] = "-x"
        elif "--orf" in self.params["redir_params"]:
            self.params["orf_to_use"] = "--orf"
        else:
            self.params["orf_to_use"] = ""
        
        
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Initialize a dict fot CARD_RGI results
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash            
            if "CARD_RGI.json" in self.sample_data[sample].keys():
                self.write_warning("That's strange, you already have a CARD_RGI slot defined for sample %s. Make sure you are not overwriting!\n" ,sample)
        
        
        sample_has_nucl = project_has_nucl = False
        if "scope" not in self.params:
            # If all samples have fasta.nucl:
            if all(map(lambda x: "fasta.nucl" in self.sample_data[x], self.sample_data["samples"])):
                sample_has_nucl = True
            if "fasta.nucl" in self.sample_data:
                project_has_nucl = True
            if sample_has_nucl and project_has_nucl:
                raise AssertionExcept("Both sample and project fasta exists. You must specify 'scope'")
            elif sample_has_nucl:
                self.params["scope"] = "sample"
            elif project_has_nucl:
                self.params["scope"] = "project"
            else:
                raise AssertionExcept("No fasta exists in either samples or project!")
            
            
        if self.params["scope"] == "sample":
            # Assert that all samples have nucleotide fasta files:
            for sample in self.sample_data["samples"]:    
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a fasta file\n", sample)
        elif self.params["scope"] == "project":
            try:
                self.sample_data["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("Project does not have a fasta file\n")
                
        
    def create_spec_wrapping_up_script(self): #not use yet
        """ Add stuff to check and agglomerate the output data
        """

        if self.params["scope"] == "sample":
            self.make_sample_file_index()   # see definition below
            
                    
            try:
                self.params["merge_script_path"]
            except KeyError:
                self.write_warning("You did not supply a merge_script_path parameter. Will not merge reports...\n")
                self.script = ""
            else:
                self.script = "%s \\\n\t--ind %s \\\n\t--output %s\n\n" % (self.params["merge_script_path"], \
                                                                            self.sample_data["CARD_RGI.files_index"], \
                                                                            self.base_dir + "CARD_RGI.merged_reports.tsv")
                self.sample_data["CARD_RGI.merged_reports"] = self.base_dir + "CARD_RGI.merged_reports.tsv"
                self.stamp_file(self.sample_data["CARD_RGI.merged_reports"])
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        if self.params["scope"] == "sample":

            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""
                
                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
                
                 
                self.script += "cd %s\n\n" % use_dir
         
                self.script += self.get_script_const()
                self.script += "-i %s \\\n\t" % self.sample_data[sample]["fasta.nucl"]  # Possibly let user decide if to use protein sequences?
                self.script += "-o %s \n\n" % (sample + ".CARD_RGI")

                # Add script to convert from JSON to tsv:
                if self.params["JSON2tsv_script"]:
                    self.script += "%s \\\n\t" % self.params["JSON2tsv_script"]
                    if self.params["orf_to_use"]:
                        self.script += "%s %s \\\n\t" % (self.params["orf_to_use"], self.params["redir_params"][self.params["orf_to_use"]])
                    self.script += "-i %s \\\n\t" % (sample + ".CARD_RGI.json")
                    self.script += "-o %s \n\n" % (sample + ".CARD_RGI")
                
                
                
                self.script += "cd %s\n\n" % self.pipe_data["home_dir"]
         

                self.sample_data[sample]["CARD_RGI.json"] = sample_dir + sample + ".CARD_RGI.json"
                self.stamp_file(self.sample_data[sample]["CARD_RGI.json"])
                
                if self.params["JSON2tsv_script"]:
                    self.sample_data[sample]["CARD_RGI.tsv"] = sample_dir + sample + ".CARD_RGI.txt"
                    self.stamp_file(self.sample_data[sample]["CARD_RGI.tsv"])

                
                

                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                        
        else:  # scope = project

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
            
             
            self.script += "cd %s\n\n" % use_dir
     
            self.script += self.get_script_const()
            self.script += "-i %s \\\n\t" % self.sample_data["fasta.nucl"]  # Possibly let user decide if to use protein sequences?
            self.script += "-o %s \n\n" % (self.sample_data["Title"] + ".CARD_RGI")

            # Add script to convert from JSON to tsv:
            if self.params["JSON2tsv_script"]:
                self.script += "%s \\\n\t" % self.params["JSON2tsv_script"]
                if self.params["orf_to_use"]:
                    self.script += "%s %s \\\n\t" % (self.params["orf_to_use"], self.params["redir_params"][self.params["orf_to_use"]])
                self.script += "-i %s \\\n\t" % (self.sample_data["Title"] + ".CARD_RGI.json")
                self.script += "-o %s \n\n" % (self.sample_data["Title"] + ".CARD_RGI")
            
            
            
            self.script += "cd %s\n\n" % self.pipe_data["home_dir"]
     

            self.sample_data["CARD_RGI.json"] = self.base_dir + self.sample_data["Title"] + ".CARD_RGI.json"
            self.stamp_file(self.sample_data["CARD_RGI.json"])
            
            if self.params["JSON2tsv_script"]:
                self.sample_data["CARD_RGI.tsv"] = self.base_dir + self.sample_data["Title"] + ".CARD_RGI.txt"
                self.stamp_file(self.sample_data["CARD_RGI.tsv"])

            
            

            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                  
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "CARD_RGI.files_index.txt", "w") as index_fh:
            index_fh.write("#Sample\tBLAST_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["CARD_RGI.tsv"]))
                
        self.sample_data["CARD_RGI.files_index"] = self.base_dir + "CARD_RGI.files_index.txt"
        
  