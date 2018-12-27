#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

""" A class defining a pipeline.

This class takes input files: samples and parameters, and creates a qsub pipeline, including dependencies
Actual work is done by calling other class types: PLCStep and PLCName
"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_GenotypeGVCFs(Step):
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


        use_dir = self.local_start(self.base_dir)
        
        for chr in self.params["chrom_list"].split(','):
            chr = chr.strip()
            my_variant_string = ""

            for cohort_gvcf in self.sample_data["cohorts"]:
                my_variant_string += "    --variant " + self.sample_data[cohort_gvcf][chr]["g.vcf"] + " \\\n"     # Getting list of samples within cohort

#            self.sample_data["original_samples"] = self.sample_data["samples"]
#            self.sample_data["samples"] = self.params["chrom_list"]

                # Name of specific script:
            
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,chr])
            self.script = ""
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(chr)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            my_string = """
                
                
                
echo '\\n---------- Create gvcf file -------------\\n'
%(GATK_path)s \\
    -T GenotypeGVCFs \\
    -R %(genome_reference)s \\
%(my_variant)s\t-L %(my_chrom)s \\
    -o %(output_vcf_creation)s
                    
            """ % { 
                    "GATK_path" : self.params["script_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "output_vcf_creation" : "{d}{s}.joint_genotyping.vcf".format(d = use_dir, s = chr), #dir????????????
                    "my_variant" : my_variant_string,
                    "my_chrom" : chr
            }
            
            
            self.script += my_string
            #self.get_script_env_path()
            
            self.sample_data[chr] = dict()
            
            self.sample_data[chr]["vcf"]= "{d}{s}.joint_genotyping.vcf".format(d = use_dir, s = chr)
            self.stamp_file(self.sample_data[chr]["vcf"])
            
            self.local_finish(use_dir,sample_dir)

            self.create_low_level_script()
