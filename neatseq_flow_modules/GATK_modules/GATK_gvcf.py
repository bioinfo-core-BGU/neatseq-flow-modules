#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

""" A class defining a pipeline.

This class takes input files: samples and parameters, and creates a qsub pipeline, including dependencies
Actual work is done by calling other class types: PLCStep and PLCName
"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_GATK_gvcf(Step):
    """ A class that defines a pipeline step name (=instance).
    """


    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

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
        
        # Prepare a list to store the qsub names of this steps scripts (will then be put in pipe_data and returned somehow)
        self.qsub_names=[]
        
       
        # Each iteration must define the following class variables:
            # spec_qsub_name
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            for chr in self.params["chrom_list"].split(','):      # Getting list of cohorts out of cohort_hash
                chr = chr.strip()
                # Name of specific script:
            
                self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample,chr])
                self.script = ""
            
            
                my_pre_string = """
cd %(sample_dir)s


echo '\\n---------- Create gvcf file -------------\\n'
"""% { "sample_dir": sample_dir }

                my_string = """-T HaplotypeCaller \\
    -R %(genome_reference)s \\
    -I %(output_duplicates)s \\
    --emitRefConfidence GVCF \\
    --variant_index_type LINEAR \\
    --variant_index_parameter 128000 \\
    -o %(output_gvcf_creation)s \\
    -L %(my_chrom)s
                        
                """ % {
                        "genome_reference" : self.params["genome_reference"],
                        "output_duplicates" : self.sample_data[sample]["bam"],
                        "output_gvcf_creation" : sample_dir + sample + "_chr_" + chr + ".g.vcf",
                        "my_chrom" : chr
                }          
                
                self.script = my_pre_string
                self.script += self.get_script_const()
                self.script += my_string

                #self.get_script_env_path()
            
                self.sample_data[sample][chr] = {}
                self.sample_data[sample][chr]["GATK_g.vcf"] = sample_dir + sample + "_chr_" + chr + ".g.vcf"

                self.stamp_file(self.sample_data[sample][chr]["GATK_g.vcf"])
            
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

                self.create_low_level_script()
                    
