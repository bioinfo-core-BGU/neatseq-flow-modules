# -*- coding: UTF-8 -*-
""" 
``add_trinity_tags`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for adding the tags required by Trinity to the ends of the read names.
See the `Strand specific assembly`_ section of the Trinity manual.

The module uses awk, so you don't need to pass a ``script_path``. 
Since you must pass a ``script_path``, just leave it blank. 


.. _Strand specific assembly: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#strand_specific_assembly
 
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    * ``fastq`` files in at least one of the following slots:
        
        * ``sample_data[<sample>]["fastq.F"]``
        * ``sample_data[<sample>]["fastq.R"]``
        * ``sample_data[<sample>]["fastq.S"]``

    
Output:
~~~~~~~~~~~~~

    * puts ``fastq`` output files (with added tags) in the following slots:
        
        * ``sample_data[<sample>]["fastq.F"]``
        * ``sample_data[<sample>]["fastq.R"]``
        * ``sample_data[<sample>]["fastq.S"]``

                
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trintags:
        module:      add_trinity_tags
        base:        trim1
        script_path: NOT_USED

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from  modules.global_defs import ZIPPED_EXTENSIONS, ARCHIVE_EXTENSIONS, KNOWN_FILE_EXTENSIONS



__author__ = "Menachem Sklarz"
__version__ = "1.1.0"



# A dict for conversion of types of sample data to positions in fasta structure:
fasta_types_dict = {"Nucleotide":"fasta.nucl","Protein":"fasta.prot"}

class Step_add_trinity_tags(Step):

    
    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = "trin_tags.fq"
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)
         

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
    def build_scripts(self):
        

        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # Adding tags to Forward and Reverse files only
            for direction in ["Forward","Reverse","Single"]:
                file_slot = "fastq." + direction[0]  # file_slot is "fastq.F", "fastq.R" and "readS" for "Forward", "Reverse" and "Single" resepctively
                if (file_slot in self.sample_data[sample].keys()):
                    self.script = ""
                    direction_tag = direction[0] # Get first letter in direction
                    # Name of specific script:
                    self.spec_script_name = "_".join([self.step,self.name,sample,direction_tag]) 
                    
                    # This line should be left before every new script. It sees to local issues.
                    # Use the dir it returns as the base_dir for this step.
                    use_dir = self.local_start(self.base_dir)

                    
                    baseFN = os.path.basename(self.sample_data[sample][file_slot])
                    # TODO: Remove ".fq" in middle of file name
                    # Setting filenames before adding output arguments to script
                    fq_fn = ".".join([baseFN, self.file_tag])  #The filename containing the end result. Used both 
                    
                    
                    
                    # self.script += self.get_script_const()
                    ""
                    self.script += "awk '{ if (NR%%4==1) { gsub(\" \",\"_\"); print $0\"%(tag)s\" } else { print } }' \\\n\t" % {"tag" : {"R":"/2","F":"/1","S":""}[direction[0]]}
                    self.script += "%s \\\n\t" % self.sample_data[sample][file_slot]
                    self.script += "> %s\n\n" % (use_dir + fq_fn)

                    
                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                    
                    # Store file in active file for sample:
                    self.sample_data[sample][file_slot] = (self.base_dir + fq_fn)
                    self.stamp_file(self.sample_data[sample][file_slot])
                    
                    self.create_low_level_script()
        