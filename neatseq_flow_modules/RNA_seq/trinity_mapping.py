# -*- coding: UTF-8 -*-
""" 
``trinity_mapping``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for running ``align_and_estimate_abundance.pl`` on a Trinity assembly and the raw reads.

Tested on versions 2.4.0 and 2.5.0 of Trinity.

See the `align_and_estimate_abundance.pl`_ script documentation.

.. _align_and_estimate_abundance.pl: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#estimating-transcript-abundance


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

* Puts output files in the following slots:
        
    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["unsorted_bam"]``  (If ``--coordsort_bam`` is passed in redirects)
    * ``sample_data[<sample>]["isoforms.results"]``
    * ``sample_data[<sample>]["genes.results"]``

Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide fasta slot should be used"
    "redirects: --gene_trans_map", "path or empty", "If empty, use internal gene_trans_map. If path, use path as gene_trans_map for all samples. If not passed, performs analysis on isoform level only"
    "redirects: --trinity_mode", "", "If set, will create a gene_trans_map for each sample and store it as sample gene_trans_map"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trin_map1:
        module:               trinity_mapping
        base:                 trinity1
        script_path:          {Vars.paths.align_and_estimate_abundance}
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
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_trinity_mapping(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "trin_mapping"
        

        if "--est_method" not in self.params["redir_params"]:
            raise AssertionExcept("You must pass an --est_method to trin_mapping.\n")

        # Is used below... 
        self.est_method = self.params["redir_params"]["--est_method"]
        if self.est_method == "kallisto":
            raise AssertionExcept("Method 'kallisto' is not defined yet!")
            # To define, find out what the per isoform and per gene output files are named and fill the names in the dictionary called file_suffix_ind, below.
        
        if "--aln_method" in self.params["redir_params"]:
            # raise AssertionExcept("You must pass an --aln_method to trin_mapping\n")
            self.params["aln_method"] = self.params["redir_params"]["--aln_method"]
            del self.params["redir_params"]["--aln_method"]
        else:
            self.params["aln_method"] = None

            
        if not self.params["aln_method"] and self.est_method.lower() in ["rsem","express"]: 
            raise AssertionExcept("For RSEM and eXpress, you must supply an 'aln_method' parameter")

        if "scope" not in self.params:
            raise AssertionExcept("Please specify a 'scope': Either 'sample' or 'project'.")
        
            
        for redir2remove in ["--transcripts", "--output_dir", "--left", "--right", "--single", "--prep_reference"]:
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
            # print self.sample_data.keys()
            if "fasta.nucl" not in self.sample_data.keys():
                raise AssertionExcept("It seems there is no project-wide assembly.")
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'.")


        
        # If "bam" required as input method, make sure a bam exists for all samples:
        if self.params["aln_method"] == "bam":
            for sample in self.sample_data["samples"]:
                if "bam" not in self.sample_data[sample]:
                    raise AssertionExcept("It seems there is no BAM file for the sample.", sample)
        
        # Dealing with gene_trans_map:
        if self.params["scope"] == "project":
            if "--gene_trans_map" in self.params["redir_params"]:
                self.use_gene_trans_map = True
                if self.params["redir_params"]["--gene_trans_map"]:  # If value was passed
                    self.sample_data["gene_trans_map"] = self.params["redir_params"]["--gene_trans_map"]
                else:  # If passed empty, use internal:
                    if "gene_trans_map" in self.sample_data:
                        self.params["redir_params"]["--gene_trans_map"] = self.sample_data["gene_trans_map"]
                    else:
                        raise AssertionExcept("Expecting 'gene_trans_map' in project but none found.\n")

            elif "--trinity_mode" in self.params["redir_params"]:
                self.sample_data["gene_trans_map"] = "%s.gene_trans_map" % self.sample_data["fasta.nucl"]
                self.use_gene_trans_map = True
            else:
                self.use_gene_trans_map = False
        else: # sample scope
            if "--gene_trans_map" in self.params["redir_params"]:
                self.use_gene_trans_map = True
                if self.params["redir_params"]["--gene_trans_map"]:  # If value was passed
                    for sample in self.sample_data["samples"]:
                        self.sample_data[sample]["gene_trans_map"] = self.params["redir_params"]["--gene_trans_map"]
                else:  # If passed empty, use internal:
                    if "gene_trans_map" in self.sample_data:
                        self.params["redir_params"]["--gene_trans_map"] = self.sample_data[sample]["gene_trans_map"]
                    else:
                        raise AssertionExcept("Expecting 'gene_trans_map' in sample but none found.\n", sample)

            elif "--trinity_mode" in self.params["redir_params"]:
                self.sample_data[sample]["gene_trans_map"] = "%s.gene_trans_map" % self.sample_data[sample]["fasta.nucl"]
                self.use_gene_trans_map = True
            else:
                self.use_gene_trans_map = False
        
         
         
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
  
        if all([self.params["scope"] == "project",  \
                self.params["aln_method"] not in ["bam", None]]):
            
            self.script = ""

       
            # Create script and write to SCRPT
            # First: transcript preparation (with --pre_reference arg)
            self.script += self.get_script_const()
            self.script += "--aln_method %s \\\n\t"   % self.params["aln_method"]
            self.script += "--transcripts %s \\\n\t"   % self.sample_data["fasta.nucl"]
            self.script += "--prep_reference \n\n"
            

            
        else:
            pass



    def build_scripts(self):
        
        
        file_suffix_ind = {
            "rsem": {
                "isoforms": "RSEM.isoforms.results",
                "genes":    "RSEM.genes.results"},
            "salmon": {
                "isoforms": "quant.sf",
                "genes":    "quant.sf.genes"},
            "kallisto": {
                "isoforms": "",
                "genes":    ""    },
            "express": {
                "isoforms": "results.xprs",
                "genes":    "results.xprs.genes"}
        }
        
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

            # # Adding single reads to end of left (=forward) reads
            # if single != "" and forward != "":
                # forward = ",".join([forward,single])

            
            transcripts = self.sample_data[sample]["fasta.nucl"] \
                if self.params["scope"] == "sample" \
                else  self.sample_data["fasta.nucl"]
                

            if all([self.params["scope"] == "sample",   \
                    "aln_method" in self.params,        \
                    self.params["aln_method"] not in ["bam", None]]):
                self.script += "# Preperaing the reference for analysis:\n\n"
                self.script += self.get_script_const()
                self.script += "--aln_method %s \\\n\t"   % self.params["aln_method"]
                self.script += "--transcripts %s \\\n\t"   % transcripts
                self.script += "--prep_reference \n\n"


                
        
            # Create script and write to SCRPT
            # First: transcript preparation (with --pre_reference arg) 
            #  - This is done with preliminary script (see create_spec_preliminary_script())
            self.script += self.get_script_const()
            
            if self.params["aln_method"] == "bam":
                self.script += "--aln_method %s \\\n\t"   % self.sample_data[sample]["bam"]  # Checked above. BAM must exist.
            elif self.params["aln_method"] == None:
                pass
            else:
                self.script += "--aln_method %s \\\n\t"   % self.params["aln_method"]

                
            self.script += "--transcripts %s \\\n\t"   % transcripts
            self.script += "--output_dir %s \\\n\t"    % use_dir
            if (forward): 
                self.script += "--left %s \\\n\t"      % forward
                self.script += "--right %s \\\n\t"       % reverse
            elif (single):
                self.script += "--single %s \\\n\t"      % single
            else:
                pass # No reads. This should be caught above...
                
            self.script = self.script.rstrip("\\\n\t") + "\n\n"

            # Stroing files:
            mv_data = {"dir" : use_dir, 
                       "src" : file_suffix_ind[self.est_method.lower()]["isoforms"],
                       "trg" : ".".join([sample,file_suffix_ind[self.est_method.lower()]["isoforms"]])}

            self.script += "mv {dir}{src} {dir}{trg}\n".format(**mv_data)
            self.sample_data[sample]["isoforms.results"] = "{dir}{trg}".format(**mv_data)
            self.stamp_file(self.sample_data[sample]["isoforms.results"])

            if self.use_gene_trans_map:  # Produce gene files:
                mv_data["src"] = file_suffix_ind[self.est_method.lower()]["genes"]
                mv_data["trg"] = ".".join([sample,file_suffix_ind[self.est_method.lower()]["genes"]])
                self.script += "mv {dir}{src} {dir}{trg}\n".format(**mv_data)
                self.sample_data[sample]["genes.results"] = "{dir}{trg}".format(**mv_data)
                self.stamp_file(self.sample_data[sample]["genes.results"])

            
            # Store bam files
            if self.est_method.lower() in ["rsem","express"]: 
                self.sample_data[sample]["bam"] = "{dir}{method}.bam".format(dir    = sample_dir, \
                                                                             method = self.params["aln_method"])
                self.stamp_file(self.sample_data[sample]["bam"])

                self.sample_data[sample]["mapper"] = "%s" % self.params["aln_method"]
                if "--coordsort_bam" in self.params["redir_params"]:
                    self.sample_data[sample]["unsorted_bam"] = self.sample_data[sample]["bam"]
                    self.stamp_file(self.sample_data[sample]["unsorted_bam"])
                    self.sample_data[sample]["bam"] = "{dir}{method}.csorted.bam".format(dir    = sample_dir, \
                                                                             method = self.params["aln_method"])
                    self.stamp_file(self.sample_data[sample]["bam"])


                
            
            self.sample_data[sample]["reference"] = transcripts
            
            
           
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
         
            self.create_low_level_script()
                        
                
            
     