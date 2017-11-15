# -*- coding: UTF-8 -*-
""" 
``htseq_count``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running htseq-count:

See htseq-count documentation `here <http://htseq.readthedocs.io/en/master/count.html>`_.



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fastq files in one of the following slots:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["sam"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts the output file in:
    ``self.sample_data[<sample>]["HTSeq.counts"]``
`

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "gff", "path to bowtie1 index", "If not given, will look for a project bowtie1 index and then for a sample bowtie1 index"
    "-f|--format", "sam | bam", "In redirects. Tells htseq-count which file to use. If not specified, will use whichever file exists."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

::

    htseq_c1:
        module:         htseq_count
        base:           samtools_STAR1
        script_path:    /storage16/app/bioinfo/python_packages/bin/htseq-count
        gtf:            /fastspace/bioinfo_databases/STAR_GRCh38_Gencode21/gencode.v21.annotation.gtf
        redirects:
            --format:   bam
            -s:         no
            -m:         intersection-nonempty



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Anders, S., Pyl, P.T. and Huber, W., 2015. **HTSeq—a Python framework to work with high-throughput sequencing data**. *Bioinformatics*, 31(2), pp.166-169.

"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_htseq_count(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "HTSeqCount"

        if not "gff" in self.params.keys():
            raise AssertionExcept("You must pass a 'gff' parameter.")

        if "-f" in self.params["redir_params"] and "--format" in self.params["redir_params"]:
            raise AssertionExcept("Please do not define both -f and --format")
        elif "--format" in self.params["redir_params"]:
            pass
        elif "-f" in self.params["redir_params"]:
            self.params["redir_params"]["--format"] = self.params["redir_params"]["-f"]
            del self.params["redir_params"]["-f"]
        else:
            self.write_warning("No format specified. Guessing...")
        
        if "--format" in self.params["redir_params"] and self.params["redir_params"]["--format"] not in ["sam","bam"]:
            raise AssertionExcept("-f/--format can be either 'bam' or 'sam'")
            

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

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
        
            
    
    
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            
                
            if "--format" in self.params["redir_params"]:
                if self.params["redir_params"]["--format"] == "sam":
                    input_file = self.sample_data[sample]["sam"]
                elif self.params["redir_params"]["--format"] == "bam":
                    input_file = self.sample_data[sample]["bam"]
                else:
                    pass
            else:
                if "sam" in self.sample_data[sample]:
                    input_file = self.sample_data[sample]["sam"]
                elif "bam" in self.sample_data[sample]:
                    input_file = self.sample_data[sample]["bam"]
                else:
                    raise AssertionExcept("No 'sam' or 'bam' exist for sample", sample)

 
            # Define location and prefix for output files:
            output_file = sample + "_htseq_count"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            
            self.script += "%s \\\n\t" % input_file
            self.script += "%s \\\n\t" % self.params["gff"]
            self.script += " > %s \n\n" % (use_dir + output_file)
            
            #save the htseq count results	
            

            self.sample_data[sample]["HTSeq.counts"] = (sample_dir + output_file)
            self.stamp_file(self.sample_data[sample]["HTSeq.counts"])
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
        
