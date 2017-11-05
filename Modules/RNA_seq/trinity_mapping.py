# -*- coding: UTF-8 -*-
""" 
``trinity_mapping``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for running ``align_and_estimate_abundance.pl`` on a Trinity assembly and the raw reads.

See the script documentation `here <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#estimating-transcript-abundance>`_.



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    * ``fastq`` files in at least one of the following slots:
        
        * ``sample_data[<sample>]["fastq.F"]``
        * ``sample_data[<sample>]["fastq.R"]``
        * ``sample_data[<sample>]["fastq.S"]``
    
    * A Trinity assembly in one of (depending on ``scope``)

        * ``sample_data[<sample>]["fasta.nucl"]``
        * ``sample_data["fasta.nucl"]``    
    
Output:
~~~~~~~~~~~~~

* Puts BAM output files in the following slots:
        
    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["isoforms.results"]``
    * ``sample_data[<sample>]["genes.results"]``

Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide fasta slot should be used"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trin_map1:
        module:               trinity_mapping
        base:                 trinity1
        script_path:          /fastspace/bioinfo_apps/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl
        redirects:
            --est_method:     RSEM
            --aln_method:     bowtie
            --trinity_mode:
            --seqType:        fq
        




References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011. **Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data**. *Nature biotechnology*, 29(7), p.644.

"""



import os
import sys
import re
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_trinity_mapping(Step):

    
    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = "trin_mapping"
        
        try:
            self.params["redir_params"]["--aln_method"]
        except KeyError:
            raise AssertionExcept("You must pass an --aln_method to trin_mapping\n")
        
        if "scope" not in self.params:
            raise AssertionExcept("Please specify a 'scope': Either 'sample' or 'project'.")
        
        
        for redir2remove in ["--transcripts", "--output_dir", "--output_prefix", "--left", "--right", "--single", "--prep_reference"]:
            if redir2remove in self.params["redir_params"]:
                del self.params["redir_params"][redir2remove]
                self.write_warning("You are not supposed to specify %s in redirects. We set it automatically" % redir2remove)
                
                
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        if self.params["scope"] == "sample":
            # Check that "fasta" and "assembly" exist (signs that trinity has been executed)
            for sample in self.sample_data["samples"]:
                if "fasta.nucl" not in self.sample_data[sample]:
                    raise AssertionExcept("It seems there is no sample-wide assembly.", sample)
        elif self.params["scope"] == "project":
            print self.sample_data.keys()
            if "fasta.nucl" not in self.sample_data.keys():
                raise AssertionExcept("It seems there is no project-wide assembly.")
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'.")

        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
  
        if self.params["scope"] == "project":
            
            self.script = ""

       
            # Create script and write to SCRPT
            # First: transcript preparation (with --pre_reference arg)
            self.script += self.get_script_const()

            self.script += "--transcripts %s \\\n\t"   % self.sample_data["fasta.nucl"]
            self.script += "--prep_reference \n\n"
        else:
            self.script = ""

            
         

    def build_scripts(self):
        
        
        # Loop over samples and concatenate read files to $forward and $reverse respectively
        # add check if paired or single !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
             
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])

            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            

            # Procedure for preparing 
            # Repeating procedure as done in trinity step:
            # If both F and R reads exist, adding them to forward and reverse
            # Assuming upstream input testing to check that if there are F reads then there are also R reads.
            # Setting variables to empty strings
            single = forward = reverse = ""

            if "fastq.F" in self.sample_data[sample]:
                forward = self.sample_data[sample]["fastq.F"]
                reverse = self.sample_data[sample]["fastq.R"]
            if "fastq.S" in self.sample_data[sample]:
                single = self.sample_data[sample]["fastq.S"]

            # Adding single reads to end of left (=forward) reads
            if single != "" and forward != "":
                forward = ",".join([forward,single])

            
            transcripts = self.sample_data[sample]["fasta.nucl"] \
                if self.params["scope"] == "sample" \
                else  self.sample_data["fasta.nucl"]
                
                
            if self.params["scope"] == "sample":  # For sample scope assembly, add prep_reference per sample.
                self.script += "# Preperaing the reference for analysis:\n\n"
                self.script += self.get_script_const()

                self.script += "--transcripts %s \\\n\t"   % transcripts
                self.script += "--prep_reference \n\n"

        
            # Create script and write to SCRPT
            # First: transcript preparation (with --pre_reference arg) 
            #  - This is done with preliminary script (see create_spec_preliminary_script())
            self.script += self.get_script_const()

            self.script += "--transcripts %s \\\n\t"   % transcripts
            self.script += "--output_dir %s \\\n\t"    % use_dir
            self.script += "--output_prefix %s \\\n\t" % sample
            if (forward): 
                self.script += "--left %s \\\n\t"      % forward
                self.script += "--right %s \n\n"       % reverse
            elif (single):
                self.script += "--single %s \n\n"      % single
 

            self.sample_data[sample]["bam"] = "%s%s.bowtie.bam" % (sample_dir,sample)
            
            # It is possible to pass a bam file instead of alignment method. 
            # This option is not covered here yet...
            self.sample_data[sample]["mapper"] = "trinity/%s" % self.params["redir_params"]["--aln_method"]
            
            self.sample_data[sample]["reference"] = transcripts
            
            self.sample_data[sample]["genes.results"] = "%s%s.genes.results" % (sample_dir,sample)
            self.sample_data[sample]["isoforms.results"] = "%s%s.isoforms.results" % (sample_dir,sample)
            
            self.stamp_file(self.sample_data[sample]["genes.results"])
            self.stamp_file(self.sample_data[sample]["isoforms.results"])
            self.stamp_file(self.sample_data[sample]["bam"])

           
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
         
                
            
            
            self.create_low_level_script()
                        
                
            
     