# -*- coding: UTF-8 -*-
"""
``merge`` :sup:`*`
-------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging <and unzipping> fastqc and fasta files

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A list of fastq files in the following slots:

    * ``sample_data[<sample>]["Forward"]``
    * ``sample_data[<sample>]["Reverse"]``
    * ``sample_data[<sample>]["Single"]``

* or a list of fasta files in the following slots:

    * ``sample_data[<sample>]["Nucleotide"]``
    * ``sample_data[<sample>]["Protein"]``
    
* or a list of BAM/SAM files in the following slots:

    * ``sample_data[<sample>]["SAM"]``
    * ``sample_data[<sample>]["BAM"]``
    * ``sample_data[<sample>]["REFERENCE"]`` - The reference fasta used to align the reads in the BAM/SAM files.
    
* or a list of VCF files in the following slots:

    * ``sample_data[<sample>]["VCF"]``
    * ``sample_data[<sample>]["G.VCF"]``
    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* puts fastq output files in the following slots:

    * ``sample_data[<sample>]["fastq.F"|"fastq.R"|"fastq.S"]``
        
* puts fasta output files in the following slots:
    
    * ``sample_data[<sample>]["fasta.nucl"|"fasta.prot"]``

* puts SAM/BAM output files in the following slots:
    
    * ``sample_data[<sample>]["sam|bam|reference"]``

* puts VCF and G.VCF output files in the following slots:
    
    * ``sample_data[<sample>]["vcf|g.vcf"]``

.. note:: In the *merge* parameters, set the *script_path* parameter according to the type of raw files you've got. 
    e.g., if they are gzipped, it should be ``gzip -cd``, etc.

.. attention:: If you want to do something more complex with the combined files, you can use the ``pipe`` parameter to send extra commands to be piped on the files after the main command. **This is an experimental feature and should be used with care**.

    e.g.: You can get files from a remote location by setting ``script_path`` to ``curl`` and ``pipe`` to ``gzip -cd``. This will download the files with curl, unzip them and concatenate them into the target file.  In the sample file, specify remote URLs instead of local pathes. **This will work only for one file per sample**.
    
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "pipe", "", "Additional commands to be piped on the files before writing to file."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    merge1:
        module: merge
        script_path: gzip -cd

::

    merge1:
        module: merge
        script_path: curl
        pipe:  gzip -cd
        
"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from  neatseq_flow.modules.global_defs import ZIPPED_EXTENSIONS, ARCHIVE_EXTENSIONS, KNOWN_FILE_EXTENSIONS

import yaml



__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_merge(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "merge"
        
        # Load YAML of file type stored in merge_file_types.yml
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"merge_file_types.yml"),"r") as fileh:
            try:
                self.default_src_trg_map = yaml.load("".join(fileh.readlines()),  Loader=yaml.SafeLoader)
            except yaml.YAMLError, exc:
                if hasattr(exc, 'problem_mark'):
                    mark = exc.problem_mark
                    print "Error position: (%s:%s)" % (mark.line+1, mark.column+1)
                    print mark.get_snippet()
                raise AssertionExcept("Error loading file types index 'merge_file_types.yml'")
            except:
                raise AssertionExcept("Error loading file types index 'merge_file_types.yml'")
            
    
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        src = set()
        # Get list of existing file types in samples file:
        for sample in self.sample_data["samples"]:
            src = src | set(self.sample_data[sample].keys())
        if "src" not in self.params:
            src = list(src)
            if "type" in src: 
                src.remove("type")   # 'type' is the type of sample (PE, etc.)
            self.params["src"] = src
        else: # Check that all src's exist somewhere in sample_data:
            if not set(self.params["src"]) < src:
                raise AssertionExcept("The following types in 'src' do not exist in samples file: %s. " % " ".join(list(set(self.params["src"]) - src)))
            
        # If script_path is not a list or is a list of length 1, extend it to length of src
        if not isinstance(self.params["script_path"],list) or len(self.params["script_path"])==1:
            self.params["script_path"] = [self.params["script_path"] for elem in self.params["src"]]
        
        # If pipe is defined, is not a list or is a list of length 1, extend it to length of src
        if "pipe" in self.params:
            if not isinstance(self.params["pipe"],list) or len(self.params["pipe"])==1:
                self.params["pipe"] = [self.params["pipe"] for elem in self.params["src"]]
            else:
                if not len(self.params["pipe"])==len(self.params["src"]):
                    raise AssertionExcept("pipe list must be the same length as src list!")

        
        # Create trg list
        if "trg" not in self.params:
            if not set(self.params["src"]) < set(self.default_src_trg_map.keys()):
                raise AssertionExcept("The following types in 'src' are not recognized: %s. Please use explicit 'trg' list in merge" % " ".join(list(set(self.params["src"]) - set(self.default_src_trg_map.keys()))))
            self.params["trg"] = [self.default_src_trg_map[src][0] for src in self.params["src"]]
            
        if "ext" not in self.params:
            if not set(self.params["src"]) < set(self.default_src_trg_map.keys()):
                raise AssertionExcept("The following types in 'src' are not recognized: %s. Please use explicit 'ext' list." % " ".join(list(set(self.params["src"]) - set(self.default_src_trg_map.keys()))))
            self.params["ext"] = [self.default_src_trg_map[src][1] if len(self.default_src_trg_map[src])>1 else src for src in self.params["src"]]
            
            
            
        if not len(self.params["script_path"])==len(self.params["src"])==len(self.params["trg"])==len(self.params["ext"]):
            raise AssertionExcept("script_path, src, trg and ext lists, if defined, must all have the same lengths!")
        
        
        
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
    def build_scripts(self):
        
        
        

        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # General comment: If there is a parallel routine for each direction (forward, reverse), add this loop	
            # if  in self.sample_data[sample].keys():

            # for type_i in range(len(self.params["src"])):
            for src_type in self.params["src"]:
            
                self.script = ""
                
                # Get index of src. Will be used to extract equivalent trg, script_path and ext.
                type_i = self.params["src"].index(src_type)
                # src_type = self.params["src"][type_i]
                
                # src_type not defined for this sample. Move on.
                if src_type not in self.sample_data[sample]:
                    continue
                    
                self.spec_script_name = "_".join([self.step,self.name,sample,src_type]) 
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                
                fq_fn = ".".join([sample, src_type, self.file_tag,self.params["ext"][type_i]])          #The filename containing the end result. Used both in script and to set reads in $sample_params


                self.script += self.params["script_path"][type_i] + " \\\n\t"
                # The following line concatenates all the files in the direction separated by a " "
                self.script += " ".join(self.sample_data[sample][src_type]) 
                self.script += " \\\n\t"
                if "pipe" in self.params:
                    self.script += "| {pipe} \\\n\t".format(pipe = self.params["pipe"][type_i])
                self.script += "> %s%s \n\n"  % (use_dir, fq_fn)

                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                
                # Store file in active file for sample:
                self.sample_data[sample][self.params["trg"][type_i]] = self.base_dir + fq_fn
                
                self.stamp_file(self.sample_data[sample][self.params["trg"][type_i]])
                
                
                self.create_low_level_script()
                