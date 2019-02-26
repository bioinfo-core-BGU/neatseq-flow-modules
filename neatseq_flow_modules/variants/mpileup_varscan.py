# -*- coding: UTF-8 -*-
""" 
``mpileup_varscan``
------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for identifying variance by running mpileup and piping it's (large) output into varscan:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BAM files in the the following slots:

    * ``sample_data[<sample>]["bam"]``

* Genome reference fasta files in the the following slot (the slot should be populated by the module that created the ``bam`` file):

    * ``sample_data[<sample>]["reference"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* If ``scope`` is set to ``sample``:

    * Puts output files in:
    
        ``sample_data[<sample>]["vcf"]``
        ``sample_data[<sample>]["variants"]`` (if ``--output-vcf`` is not redirected in ``redirects``)

* If ``scope`` is set to ``project``:

    * Puts output files in:
    
        ``sample_data["vcf"]``
        ``sample_data["variants"]`` (if ``--output-vcf`` is not redirected in ``redirects``)


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "mpileup_path",  "path", "The full path to the mpileup program. You can append additional mpileup arguments after the path (see example lines)"
    "script_path",  "path", "The full path to the relevant varscan program (see example lines)."
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    mpileup_varscan1:
        module: mpileup_varscan
        base: samtools1
        script_path: /path/to/java -jar /path/to/VarScan.v2.3.9.jar mpileup2snp
        mpileup_path: /path/to/samtools mpileup --max-depth 6000
        scope: sample
        redirects:
            --min-coverage: 4
            --output-vcf:
            --variants: 1
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G. and Durbin, R., 2009. **The sequence alignment/map format and SAMtools**. *Bioinformatics*, 25(16), pp.2078-2079.

- Koboldt, D.C., Chen, K., Wylie, T., Larson, D.E., McLellan, M.D., Mardis, E.R., Weinstock, G.M., Wilson, R.K. and Ding, L., 2009. **VarScan: variant detection in massively parallel sequencing of individual and pooled samples**. *Bioinformatics*, 25(17), pp.2283-2285.
"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_mpileup_varscan(Step):
    
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        
        if "scope" not in list(self.params.keys()):
            raise AssertionExcept("You must supply a 'scope' param. Either 'sample' or 'project'")
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        for sample in self.sample_data["samples"]:
            try:
                self.sample_data[sample]["reference"]
                self.sample_data[sample]["bam"]
            except KeyError:
                raise AssertionExcept("Sample does not have a bam or reference slot in the mapping slot", sample)
                
            # If scope is sample, open a variants slot for the sample:
            if self.params["scope"] == "sample":
                try:
                    self.sample_data[sample]["variants"]
                except KeyError:
                    self.sample_data[sample]["variants"] = dict()
           
        # If scope is project, open variants slot for project
        if self.params["scope"] == "project":
            try:
                self.sample_data["project_data"]["variants"]
            except KeyError:
                self.sample_data["project_data"]["variants"] = dict()
        
            
            
        
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

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""
            
            
            #################
            ##################
            ##
            ## TODO
            ##
            ## 1. Output sample list into file in data/mpileup_varscan
            ## 2. Add parameter for user to pass mpileup location + arguments
            ## 3. Create the mpileup command 
            ## 4. pipe the output to varscan. (set script_path to location of varscan)
            ## 5. pipe varscan output to file location
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            ### Getting location of reference
            # Get list of reference fasta files from samples, and convert to set, removing duplicates
            reference_fasta = set([self.sample_data[sample]["reference"] for sample in self.sample_data["samples"]])
            # If there are more than one reference_fasta, exit. This is really really weird and should not happen
            if len(reference_fasta) > 1:
                raise AssertionExcept("There is more than one reference file for different samples. Weird!!!")
            # Convert set into list and return first, and only, element:
            reference_fasta = list(reference_fasta)[0]

            
            ### Create file with list of sample names, one per line
            with open(self.base_dir + "sample_list.txt", "w") as smp_lst:
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    smp_lst.write("%s\n" % sample)

            # Define output_suffix depending on redir_params:
            # If either --VCF or -v are passed by user
            if "--output-vcf" in self.params["redir_params"]:
                output_suffix = "vcf"
            else:
                output_suffix = "unknown"

                
            # Get constant part of script:
            self.script += "%s \\\n\t" % self.params["mpileup_path"]
            # Reference file:
            self.script += "-f %s \\\n\t" % reference_fasta
            # BAM files:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s \\\n\t" % self.sample_data[sample]["bam"]
            # Remove extra stuff from end of script:
            self.script = self.script.rstrip("\\\n\t")
            # self.script = self.script.rstrip()
            self.script += " | \\\n"
            
            self.script += "%s \\\n\t" % self.params["script_path"]
            self.script += self.get_redir_parameters_script()
            self.script += "--vcf-sample-list %s \\\n\t" % (self.base_dir+"sample_list.txt")
            self.script += "> %s\n\n"  % (use_dir + ".".join([self.sample_data["Title"],output_suffix]))
            
            
            
            
            if output_suffix == "vcf":
                self.sample_data["project_data"]["vcf"] = (self.base_dir + ".".join([self.sample_data["Title"],output_suffix]))
                self.stamp_file(self.sample_data["project_data"]["vcf"])
            else:
                self.sample_data["project_data"]["variants"] = (self.base_dir + ".".join([self.sample_data["Title"],output_suffix]))
                self.stamp_file(self.sample_data["project_data"]["variants"])


        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                                
        else:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
         
                #################
                ##################
                ##
                ## TODO
                ##
                ## 1. Output sample list into file in data/mpileup_varscan
                ## 2. Add parameter for user to pass mpileup location + arguments
                ## 3. Create the mpileup command 
                ## 4. pipe the output to varscan. (set script_path to location of varscan)
                ## 5. pipe varscan output to file location
                
                
                # Get list of reference fasta files from samples, and convert to set, removing duplicates
                reference_fasta = self.sample_data[sample]["reference"] 
                
                
                ### Create file with list of sample names, one per line
                with open(use_dir + "sample_list.txt", "w") as smp_lst:
                    smp_lst.write("%s\n" % sample)

                # Define output_suffix depending on redir_params:
                # If either --VCF or -v are passed by user
                if "--output-vcf" in self.params["redir_params"]:
                    output_suffix = "vcf"
                else:
                    output_suffix = "unknown"

                    
                # Get constant part of script:
                self.script += "%s \\\n\t" % self.params["mpileup_path"]
                # Reference file:
                self.script += "-f %s \\\n\t" % reference_fasta
                # BAM file:
                self.script += "%s \\\n\t" % self.sample_data[sample]["bam"]
                # Remove extra stuff from end of script:
                self.script = self.script.rstrip("\\\n\t")
                # self.script = self.script.rstrip()
                self.script += " | \\\n"
                
                self.script += "%s \\\n\t" % self.params["script_path"]
                self.script += self.get_redir_parameters_script()
                self.script += "--vcf-sample-list %s \\\n\t" % (use_dir + "sample_list.txt")
                self.script += "> %s%s_%s.%s\n\n"  % (use_dir, sample, self.get_step_name(),output_suffix)
                
                
                
                
                if output_suffix == "vcf":

                    self.sample_data[sample]["vcf"] = "%s%s_%s.vcf" % (use_dir, sample, self.get_step_name())
                    self.sample_data[sample]["vcf.source"] = "varscan"
                    self.stamp_file(self.sample_data[sample]["vcf"])
                            
                else:
                    self.sample_data[sample]["variants"] = "%s%s_%s.%s" % (use_dir, sample, self.get_step_name(),output_suffix)
                    self.sample_data[sample]["variants.source"] = "varscan"
                    self.stamp_file(self.sample_data[sample]["variants"])


            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                    