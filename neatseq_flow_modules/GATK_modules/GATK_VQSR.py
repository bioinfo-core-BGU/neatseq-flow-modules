#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

""" A class defining a pipeline.

This class takes input files: samples and parameters, and creates a qsub pipeline, including dependencies
Actual work is done by calling other class types: PLCStep and PLCName
"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_GATK_VQSR(Step):
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
######################################################## SNP

        use_dir = self.local_start(self.base_dir)
        for chr in self.params["chrom_list"].split(','):
            chr = chr.strip()
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],chr])

            input_file = self.sample_data[chr]["vcf"]
            recalFile_SNP = use_dir + "recalibrate_SNP_" + chr + ".recal"
            tranchesFile_SNP = use_dir + "recalibrate_SNP_" + chr + ".tranches"
            rscriptFile_SNP = use_dir + "recalibrate_SNP_plots_" + chr + ".R"
            output_AppllyRecal_SNP = use_dir + "recalibrated_snps_raw_indels_" + chr + ".vcf"
            recalFile_INDEL = use_dir + "recalibrate_INDEL_" + chr + ".recal"
            tranchesFile_INDEL = use_dir + "recalibrate_INDEL_" + chr + ".tranches"
            rscriptFile_INDEL = use_dir + "recalibrate_INDEL_plots_" + chr + ".R"
            output_AppllyRecal_INDEL = use_dir + "recalibrated_variants_" + chr + ".vcf"
            ts_filter_level_SNP = self.params["ts_filter_level_SNP"]
            ts_filter_level_INDEL = self.params["ts_filter_level_INDEL"]
            my_SNP_VariantRecalibrator_string = """
                
                
                
    echo '\\n---------- Create SNP VariantRecalibrator file -------------\\n'
    %(GATK_path)s \\
        -T VariantRecalibrator \\
        -R %(genome_reference)s \\
        -mode SNP \\
        -an QD \\
        -an FS \\
        -an SOR \\
        -an MQ \\
        -an MQRankSum \\
        -an ReadPosRankSum \\
        -an InbreedingCoeff \\
        -input %(input_file)s \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        -recalFile %(recalFile)s \\
        -tranchesFile %(tranchesFile)s \\
        -rscriptFile %(rscriptFile)s \\
    """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"],
                        "input_file" : input_file,
                        "recalFile" : recalFile_SNP,
                        "tranchesFile" : tranchesFile_SNP,
                        "rscriptFile" : rscriptFile_SNP
                }      
            my_SNP_ApplyRecalibration_string = """
                    
                    
                    
    echo '\\n---------- Create SNP ApplyRecalibration file -------------\\n'
    %(GATK_path)s \\
        -T ApplyRecalibration \\
        -R %(genome_reference)s \\
        -mode SNP \\
        -input %(input_file)s \\
        --ts_filter_level %(ts_filter_level_SNP)s \\
        -recalFile %(recalFile)s \\
        -tranchesFile %(tranchesFile)s \\
        -o %(output_AppllyRecal_SNP)s
    """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"],
                        "input_file" : input_file,
                        "ts_filter_level_SNP" : ts_filter_level_SNP,
                        "recalFile" : recalFile_SNP,
                        "tranchesFile" : tranchesFile_SNP,
                        "output_AppllyRecal_SNP" : output_AppllyRecal_SNP
                }             

    ######################################################################### INDEL
            my_INDEL_VariantRecalibrator_string = """
                    
                    
                    
    echo '\\n---------- Create INDEL VariantRecalibrator file -------------\\n'
    %(GATK_path)s \\
        -T VariantRecalibrator \\
        -R %(genome_reference)s \\
        -mode INDEL \\
        -an QD \\
        -an FS \\
        -an SOR \\
        -an MQRankSum \\
        -an ReadPosRankSum \\
        -an InbreedingCoeff \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        --maxGaussians 4 \\
        -input %(input_file)s \\
        -recalFile %(recalFile)s \\
        -tranchesFile %(tranchesFile)s \\
        -rscriptFile %(rscriptFile)s \\
    """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"],
                        "input_file" : output_AppllyRecal_SNP,
                        "recalFile" : recalFile_INDEL,
                        "tranchesFile" : tranchesFile_INDEL,
                        "rscriptFile" : rscriptFile_INDEL
                        }
                
            my_INDEL_ApplyRecalibration_string = """
                    
                    
                    
    echo '\\n---------- Create INDEL ApplyRecalibration file -------------\\n'
    %(GATK_path)s \\
        -T ApplyRecalibration \\
        -R %(genome_reference)s \\
        -mode INDEL \\
        -input %(input_file)s \\
        --ts_filter_level %(ts_filter_level_INDEL)s \\
        -recalFile %(recalFile)s \\
        -tranchesFile %(tranchesFile)s \\
        -o %(output_AppllyRecal_INDEL)s
    """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"],
                        "input_file" : output_AppllyRecal_SNP,
                        "ts_filter_level_INDEL" : ts_filter_level_INDEL,
                        "recalFile" : recalFile_INDEL,
                        "tranchesFile" : tranchesFile_INDEL,
                        "output_AppllyRecal_INDEL" : output_AppllyRecal_INDEL
                }

            for resource in self.params["resource_SNP"]:
                my_SNP_VariantRecalibrator_string = my_SNP_VariantRecalibrator_string + "    -resource:" + resource +  " \\\n"
            my_SNP_VariantRecalibrator_string = my_SNP_VariantRecalibrator_string
            for resource in self.params["resource_INDEL"]:
                my_INDEL_VariantRecalibrator_string = my_INDEL_VariantRecalibrator_string + "    -resource:" + resource +  " \\\n"
            my_INDEL_VariantRecalibrator_string = my_INDEL_VariantRecalibrator_string
            my_per_chr_string = ""
            self.script = ""
            
             
            self.script = my_SNP_VariantRecalibrator_string + my_SNP_ApplyRecalibration_string + my_INDEL_VariantRecalibrator_string + my_INDEL_ApplyRecalibration_string


            self.sample_data[chr]["vcf"]= output_AppllyRecal_INDEL
            self.stamp_file(self.sample_data[chr]["vcf"])
            
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
                    

