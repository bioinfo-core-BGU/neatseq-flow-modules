""" 
``spades_assembl`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A class that defines a module for assembling reads using spades assembler.

 
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        * fastq files in at least one of the following slots:
        
            * ``sample_data[<sample>]["fastq.F"]``
            * ``sample_data[<sample>]["fastq.R"]``
            * ``sample_data[<sample>]["fastq.S"]``
    
    
Output:
~~~~~~~~~~~~~


    * puts fasta output files in the following slots:
        
        * for sample-wise assembly:
        
            * ``sample_data[<sample>]["fasta.nucl"]``
            * ``sample_data[<sample>]["spades_assembl.contigs"]``
            * ``sample_data[<sample>]["spades_assembl.scaffolds"]``
        
        * for mega assembly (not defined yet):
        
            * ``sample_data["fasta.nucl"]``
            * ``sample_data["spades_assembl.contigs"]``
            * ``sample_data["spades_assembl.scaffolds"]``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide fasta slot should be used"
    "truncate_names", , "truncates contig names, *e.g.* '>NODE_82_length_18610_cov_38.4999_ID_165' will be changed to '>NODE_82_length_18610'"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    spades1:
        module: spades_assembl
        base: trim1
        script_path: /path/to/bin/spades.py
        truncate_names: 
        redirects:
            --careful: 

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D. and Pyshkin, A.V., 2012. **SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing**. *Journal of computational biology*, 19(5), pp.455-477.
"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_spades_assembl(Step):

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".spades.out"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)
         
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                raise AssertionExcept("project wide scope is not defined yet\n")

            elif self.params["scope"]=="sample":
                pass
                
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
        
        
            
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
    
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        if self.params["scope"] == "project":
            
            # Not defined yet
            raise AssertionExcept("project wide scope is not defined yet\n")
        
        else:
        
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

                self.script += self.get_script_const()
                self.script += "-o %s \\\n\t" % sample_dir


                if "PE" in self.sample_data[sample]["type"]:
                    self.script += "--pe1-1 %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                    self.script += "--pe1-2 %s \n\n" % self.sample_data[sample]["fastq.R"]
                elif "SE" in self.sample_data[sample]["type"]:
                    self.script += "--s1 %s \n\n" % self.sample_data[sample]["fastq.S"]
                elif "PE" in self.sample_data[sample]["type"] and "SE" in self.sample_data[sample]["type"]:       # Mixed!!
                    self.script += "--pe1-1 %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                    self.script += "--pe1-2 %s \\\n\t" % self.sample_data[sample]["fastq.R"]
                    self.script += "--s1 %s \n\n" % self.sample_data[sample]["fastq.S"]
                else:
                    raise AssertionExcept("Strange type configuration for sample\n" ,sample)
                    
                # For prokka compliance, you can request a truncation of the contig names
                # e.g. ">NODE_82_length_18610_cov_38.4999_ID_165" will be changed to ">NODE_82_length_18610"
                if "truncate_names" in self.params.keys():
                    self.script += """
# Truncating contig names for prokka compliance
cat %(contigs)s  | cut -f 1-2 -d '_' > %(shortnames)s
mv -f %(shortnames)s %(contigs)s \n\n""" % {"contigs":sample_dir + "contigs.fasta", "shortnames":sample_dir + "contigs.shortIDs.fasta"}
                        
                # Store results to fasta and assembly slots:
                self.sample_data[sample]["fasta.nucl"]  = sample_dir + "contigs.fasta"
                self.sample_data[sample][self.get_step_step() + "_contigs"] = sample_dir + "contigs.fasta"
                self.sample_data[sample][self.get_step_step() + "_scaffolds"] = sample_dir + "scaffolds.fasta"
                self.sample_data[sample]["assembler"] = self.get_step_step()

                self.stamp_file(self.sample_data[sample][self.get_step_step() + "_scaffolds"])
                self.stamp_file(self.sample_data[sample][self.get_step_step() + "_contigs"])

                    
                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                
                self.create_low_level_script()
                        

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        pass