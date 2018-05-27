# -*- coding: UTF-8 -*-
""" 
``CollectRnaSeqMetrics`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for RNA_seq assembly annotation using `Trinotate`_.

.. _Trinotate: https://trinotate.github.io/

.. Note:: This module will be updated in the future to support uploading of other sources of information such as RNAMMER output. See `Trinotate`_ documentation.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* A transcripts file in 
    * self.sample_data["transcripts.fasta.nucl"],
* A gene to transcript mapping file in: (produced by ``Trinity_gene_to_trans_map`` module)
    * self.sample_data["gene_trans_map"], 
* A protein fasta file (produced by ``TransDecoder``)
    * self.sample_data["fasta.prot"])
* Results of ``blastx`` of protein file against swissprot database:
    * self.sample_data["blast.prot"], 
* Results of ``blastn`` of transcripts file against swissprot database:
    * self.sample_data["blast.nucl"], 
* Results of ``hmmscan`` of protein file against pfam database:
    * self.sample_data["hmmscan.prot"])
    
.. Attention:: If ``scope`` is set to ``sample``, all of the above files should be in the sample scope!

    
Output:
~~~~~~~~~~~~~

* puts Trinotate report file in:

    * ``sample_data[<sample>]["trino.rep"]`` (``scope = sample``)
    * ``sample_data["trino.rep"]`` (``scope = project``)
        

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", ""
    "sqlitedb","","Path to Trinotate sqlitedb"
    "cp_sqlitedb","","Create local copy of the sqlitedb, before loading teh data (recommended)"
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_Trinotate:
        module:             Trinotate
        base:               
                            - trino_blastp_sprot
                            - trino_blastx_sprot
                            - trino_hmmscan1
        script_path:        {Vars.paths.Trinotate}
        scope:              project
        sqlitedb:           {Vars.databases.trinotate.sqlitedb}
        cp_sqlitedb:    
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_CollectRnaSeqMetrics(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        
        self.arg_separator='='
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        if "scope" not in self.params:
            raise AssertionExcept("No 'scope' specified.")
        elif self.params["scope"]=="project":
            if "bam" not in self.sample_data and "sam" not in self.sample_data:
                raise AssertionExcept("Project does not have BAM or SAM files.")

        elif self.params["scope"]=="sample":
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                if "bam" not in self.sample_data[sample] and "sam" not in self.sample_data[sample]:
                    raise AssertionExcept("Sample does not have BAM or SAM files.", sample)
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        
         
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):
    
        if self.params["scope"] == "project":
            self.build_scripts_project()
        else:
            self.build_scripts_sample()
            
            
    def build_scripts_project(self):
        
        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        output_basename = "{title}.RNAseqMetrics.out".format(title = self.sample_data["Title"])

        if "bam" in self.sample_data:
            input = self.sample_data["bam"]
        else:
            input = self.sample_data["sam"]
            
        # Get constant part of script:
        self.script += self.get_script_const()
        for other_file in ["REF_FLAT","RIBOSOMAL_INTERVALS"]:
            if other_file not in self.params["redir_params"] and other_file in self.sample_data:
                self.script += "{filetype}={filename} \\\n\t".format(filetype=other_file, filename=self.sample_data[other_file])
        self.script += "I={input} \\\n\tO={output}\n\n".format(input=input, output=use_dir+output_basename)
                
        # Store results to fasta and assembly slots:
        self.sample_data["RNAseqMetrics"] = self.base_dir + output_basename
        self.stamp_file(self.sample_data["RNAseqMetrics"])
        for other_file in ["REF_FLAT","RIBOSOMAL_INTERVALS"]:
            if other_file in self.params["redir_params"]:
                self.sample_data[other_file] = self.params["redir_params"][other_file]

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            output_basename = "{title}.RNAseqMetrics.out".format(title = sample)

            if "bam" in self.sample_data[sample]:
                input = self.sample_data[sample]["bam"]
            else:
                input = self.sample_data[sample]["sam"]
                
            # Get constant part of script:
            self.script += self.get_script_const()
            for other_file in ["REF_FLAT","RIBOSOMAL_INTERVALS"]:
                if other_file not in self.params["redir_params"] and other_file in self.sample_data[sample]:
                    self.script += "{filetype}={filename} \\\n\t".format(filetype=other_file, filename=self.sample_data[sample][other_file])
            self.script += "I={input} \\\n\tO={output}\n\n".format(input=input, output=use_dir+output_basename)
                    
            # Store results to fasta and assembly slots:
            self.sample_data[sample]["RNAseqMetrics"] = self.base_dir + output_basename
            self.stamp_file(self.sample_data[sample]["RNAseqMetrics"])
            for other_file in ["REF_FLAT","RIBOSOMAL_INTERVALS"]:
                if other_file in self.params["redir_params"]:
                    self.sample_data[sample][other_file] = self.params["redir_params"][other_file]

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            self.create_low_level_script()
            
     