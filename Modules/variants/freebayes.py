# -*- coding: UTF-8 -*-
""" 
``freebayes``
-----------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for identifying variants by running freebayes:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BAM files in the the following slots:

    * ``sample_data[<sample>]["bam"]``

* Genome reference fasta files in the the following slot (the slot should be populated by the module that created the ``bam`` file):

    * ``sample_data[<sample>]["reference"]``
    
..Note:: Do not specify the reference (-f), since it is filled in automatically by neatseq-flow

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* If ``scope`` is set to ``sample``:

    * Puts output files in:
    
        ``sample_data[<sample>]["vcf"]``  (if ``output_type`` is set to ``vcf``)
        ``sample_data[<sample>]["gvcf"]`` (if ``output_type`` is set to ``gvcf``)

* If ``scope`` is set to ``project``:

    * Puts output files in:
    
        ``sample_data["vcf"]``  (if ``output_type`` is set to ``vcf``)
        ``sample_data["gvcf"]`` (if ``output_type`` is set to ``gvcf``)


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "output_type",  "vcf|gvcf", "The type of output produced by freebayes. (Can be specified alternatively with appropriate redirects)"
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    freebayes1:
        module: freebayes
        base: samtools1
        script_path: /path/to/freebayes
        scope: sample
        output_type: vcf
        redirects: 
            --strict-vcf:
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Marth GT, Korf I, Yandell MD, Yeh RT, Gu Z, Zakeri H, Stitziel NO, Hillier L, Kwok PY, Gish WR: **A general approach to single-nucleotide polymorphism discovery**. *Nat Genet.* 1999, 23: 452-456. 10.1038/70570.


"""


import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_freebayes(Step):
    """ A module for running freebayes:
        requires:
            sample_data[sample]["fastq"]["mapping"]["bam"]
            sample_data[sample]["fastq"]["mapping"]["reference"]
            
            output:
            puts output files in the following slots, respectively:
            self.sample_data["variants"]["VCF"]
    """
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        if "scope" not in self.params.keys():
            raise AssertionExcept("You must supply a 'scope' param. Either 'sample' or 'project'")
            
        if "output_type" in self.params:
            if self.params["output_type"] not in ["vcf","gvcf"]:
                raise AssertionExcept("'output_type' must be either 'vcf' or 'gvcf'")
            # Remove output options from redir_params
            if filter(lambda x: x in ["--vcf","-v","--gvcf"], self.params["redir_params"]):
                self.write_warning("Removing output type from redirects")
            try:
                del self.params["redir_params"]["--vcf"]
                del self.params["redir_params"]["--v"]
                del self.params["redir_params"]["--gvcf"]
            except KeyError:
                pass
        else:
            if "--vcf" in self.params["redir_params"] or "-v" in self.params["redir_params"]:
                self.params["output_type"] = "vcf"
            elif "--gvcf" in self.params["redir_params"]:
                self.params["output_type"] = "gvcf"
            else:
                raise AssertionExcept("You did not define an output type in 'output_type' nor in redirects...")
            # Removing output parameters from redirects:
            try:
                del self.params["redir_params"]["--vcf"]
                del self.params["redir_params"]["--v"]
                del self.params["redir_params"]["--gvcf"]
            except KeyError:
                pass

            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        for sample in self.sample_data["samples"]:
            try:
                self.sample_data[sample]["reference"]
                self.sample_data[sample]["bam"]
            except KeyError:
                raise AssertionExcept("Sample does not have bam or reference file types", sample)
                
           
        
            
        
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

        if self.params["scope"] == "project":
            # Name of specific script:
            self.set_spec_script_name() #"_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Define location and prefix for output files:
            # output_prefix = sample + "_bowtie2_map"

            # Get list of reference fasta files from samples, and convert to set, removing duplicates
            reference_fasta = set([self.sample_data[sample]["reference"] for sample in self.sample_data["samples"]])
            # If there are more than one reference_fasta, exit. This is really really weird and should not happen
            if len(reference_fasta) > 1:
                raise AssertionExcept("There is more than one reference file for the samples. Weird!!!" )

            # Convert set into list and return first, and only, element:
            reference_fasta = list(reference_fasta)[0]
            
            # Get constant part of script:
            self.script += self.get_script_const()
            # Reference file:
            self.script += "-f %s \\\n\t" % reference_fasta
            # BAM files:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "-b %s \\\n\t" % self.sample_data[sample]["bam"]
            
            if self.params["output_type"] == "vcf":
                self.script += "--vcf %s%s_%s.vcf \n\n" % (use_dir, self.sample_data["Title"], self.get_step_name())

                self.sample_data["vcf"] = "%s%s_%s.vcf" % (self.base_dir,self.sample_data["Title"], self.get_step_name())
                self.sample_data["vcf.source"] = "freebayes"
                self.stamp_file(self.sample_data["vcf"])
            else:   # output_type = "gvcf"
                self.script += "--gvcf %s%s_%s.gvcf \n\n" % (use_dir, self.sample_data["Title"], self.get_step_name())

                self.sample_data["gvcf"] = "%s%s_%s.gvcf" % (self.base_dir,self.sample_data["Title"], self.get_step_name())
                self.sample_data["gvcf.source"] = "freebayes"
                self.stamp_file(self.sample_data["gvcf"])
            
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                     
        else:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.set_spec_script_name(sample=True) #"_".join([self.step,self.name,sample])
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
         
                

                # Get list of reference fasta files from samples, and convert to set, removing duplicates
                reference_fasta = self.sample_data[sample]["reference"] 
                
                # Get constant part of script:
                self.script += self.get_script_const()
                # Reference file:
                self.script += "-f %s \\\n\t" % reference_fasta
                # BAM files:
                self.script += "-b %s \\\n\t" % self.sample_data[sample]["bam"]
                
                if self.params["output_type"] == "vcf":
                    self.script += "--vcf %s%s_%s.vcf \n\n" % (use_dir, sample, self.get_step_name())
                    self.sample_data[sample]["vcf"] = "%s%s_%s.vcf" % (sample_dir, sample, self.get_step_name())
                    self.sample_data[sample]["vcf.source"] = "freebayes"
                    self.stamp_file(self.sample_data[sample]["vcf"])
                else:   # output_type = "gvcf"
                    self.script += "--gvcf %s%s_%s.gvcf \n\n" % (use_dir, sample, self.get_step_name())
                    self.sample_data[sample]["gvcf"] = "%s%s_%s.gvcf" % (sample_dir, sample, self.get_step_name())
                    self.sample_data[sample]["gvcf.source"] = "freebayes"
                    self.stamp_file(self.sample_data[sample]["gvcf"])

                

            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                         
