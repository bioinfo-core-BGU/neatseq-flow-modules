# -*- coding: UTF-8 -*-
""" 
``quast`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running quast on fasta assemblies:

QUAST is executed on the fasta file along the following lines:

* If 'scope' is specified, the appropriate fasta will be used. An error will occur if the fasta does not exist.
* If 'scope' is not specified, if a project-wide fasta exists, it will be used. Otherwise, sample-wise fasta files will be used. If none exist, an error will occur.

.. Note:: With ``compare_mode``, you tell the module to run **quast** on multiple assemblies. This is done in one of three ways:

    * If ``scope`` is sample and a single base step defined, will compare between the samples.
    * If ``scope`` is sample and there is more than one base step defined, will compare between the assemblies found in the base steps for each sample separately.
    * If ``scope`` is project, will compare between the assemblies found in the base steps at the project level.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output directory in one of:
    * ``self.sample_data["quast"]``
    * ``self.sample_data[<sample>]["quast"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Indicates whether to use a project or sample contigs file."
    "compare_mode", "", "If 'scope' is 'sample', specifies wether to analyse each sample separately or to create a single comparison report for all samples."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. A quast report for each sample separately: 

::

    quast1:
        module: quast
        base: spades1
        script_path: /path/to/quast.py
        scope: sample
        redirects:
            --fast: 

2. A quast report comparing the sample assemblies: 

::

    quast1:
        module: quast
        base: spades1
        script_path: /path/to/quast.py
        compare_mode: 
        scope: sample
        redirects:
            --fast: 

3. A quast report comparing the project assemblies from different stages of the analysis: 

::

    quast1:
        module: quast
        base: 
            - spades1
            - megahit1
        script_path: /path/to/quast.py
        compare_mode: 
        scope: project
        redirects:
            --fast: 

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gurevich, A., Saveliev, V., Vyahhi, N. and Tesler, G., 2013. **QUAST: quality assessment tool for genome assemblies**. *Bioinformatics*, 29(8), pp.1072-1075.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept
from pprint import pprint as pp


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_quast(Step):

    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".quast.out"


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" in self.params.keys():
            if self.params["scope"] == "project":
                try:  # Is there a mega-assembly?
                    self.sample_data["fasta.nucl"]
                except KeyError:   # No. Check if all samples have assemblies:
                    raise AssertionExcept("No project wide assembly!")
                else:
                    pass

                if "compare_mode" in self.params.keys():
                    self.write_warning("Ignoring 'compare_mode' in project scope")
            
            elif self.params["scope"] == "sample":
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                    # Make sure each sample has a ["fasta.nucl"] slot 
                    try:
                        self.sample_data[sample]["fasta.nucl"]
                    except KeyError:
                        raise AssertionExcept("You are trying to run QUAST with no assembly.\n" , sample)
                    else:
                        pass
            else:
                raise AssertionExcept("'scope' must be either 'project' or 'sample'")
            
        
        
        
        else:
            self.write_warning("'scope' not passed. Will try guessing...")

            try:  # Is there a mega-assembly?
                self.sample_data["fasta.nucl"]
            except KeyError:   # No. Check if all samples have assemblies:
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                    # Make sure each sample has a ["fasta.nucl"] slot 
                    try:
                        self.sample_data[sample]["fasta.nucl"]

                    except KeyError:
                        raise AssertionExcept("You are trying to run QUAST with no assembly.\n" , sample)
                
                self.params["scope"] = "sample"
                
            else:
                self.write_warning("There is a project-wide assembly. Using it.\n")
        
                self.params["scope"] = "project"
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        

        # Two options: 
        # 1. If a single base, compare between samples.
        # 2. If more than one base, compare for each sample         
        
        multiple_bases = len(self.params["base"]) > 1
        
        if self.params["scope"] == "sample": # Requested for mega-assembly

            if "compare_mode" in self.params.keys() and not multiple_bases:    # Compare sample assemblies
                # print "in here: multiple_bases"
                # print multiple_bases
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
                self.script = ""

                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                    
                # Define output filename 
                # output_filename = "".join([use_dir , sample , self.file_tag])

                self.script += self.get_script_const()
                # All other parameters are redirected!
                

                self.script += "--output-dir %s \\\n\t" % use_dir
                self.script += "--labels %s \\\n\t" % ",".join(self.sample_data["samples"])
                
                # Input file:
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    self.script += "%s \\\n\t" % self.sample_data[sample]["fasta.nucl"]

                self.script = self.script.rstrip("\\\n\t")
                self.script += "\n\n"
                

            
                self.sample_data["quast"] = self.base_dir
            

                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                          
                
                self.create_low_level_script()

            else:       # Separate quast run for each sample
                # print "in here: multiple_bases:"
                # print multiple_bases
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
                    # output_filename = "".join([use_dir , sample , self.file_tag])

                    self.script += self.get_script_const()

                    # All other parameters are redirected!
                    self.script += "--output-dir %s \\\n\t" % sample_dir
                    
                    # Input file:
                    # self.script += "%s \n\n" % self.sample_data[sample]["fasta.nucl"]
                    if multiple_bases and "compare_mode" in self.params.keys():   # More than one base
                        self.script += "--labels %s \\\n\t" % ",".join(self.params["base"])
                        for base in self.params["base"]:
                            # print base
                            try:
                                self.script += "%s \\\n\t" % self.base_sample_data[base][sample]["fasta.nucl"]
                            except:
                                raise AssertionExcept("'fasta.nucl' file for sample not found in base %s" % base, sample)
                    else:
                        self.script += "%s \n\n" % self.sample_data[sample]["fasta.nucl"]
                    

                
                    # Store BLAST result file:
                    self.sample_data[sample]["quast"] = sample_dir

                    # Wrapping up function. Leave these lines at the end of every iteration:
                    self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                              
                    
                    self.create_low_level_script()

        else:    # 'scope' = project
            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""

            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
                
                
            # Define output filename 
            # output_filename = "".join([use_dir , sample , self.file_tag])

            self.script += self.get_script_const()

            # All other parameters are redirected!
            self.script += "--output-dir %s \\\n\t" % use_dir
            
            
            # Input file:
            if multiple_bases and "compare_mode" in self.params.keys():   # More than one base
                self.script += "--labels %s \\\n\t" % ",".join(self.params["base"])
                pp( self.params["base"])
                for base in self.params["base"]:
                    print base
                    pp(self.base_sample_data[base].keys())
                    try:
                        self.script += "%s \\\n\t" % self.base_sample_data[base]["fasta.nucl"]
                    except:
                        raise AssertionExcept("'fasta.nucl' file not found in base %s" % base)
            else:
                self.script += "%s \n\n" % self.sample_data["fasta.nucl"]
            

        
            # Store BLAST result file:
            self.sample_data["quast"] = self.base_dir


            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                      
            
            self.create_low_level_script()

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        pass