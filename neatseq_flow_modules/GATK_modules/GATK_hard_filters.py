#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_GATK_hard_filters(Step):
    """ A class that defines a pipeline step name (=instance).
    """


    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "csh" by inheriting instances

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
            raw_snps = use_dir + "raw_snps_" + chr + ".vcf"
            raw_indel = use_dir + "raw_indel_" + chr + ".vcf"
            filtered_snps = use_dir + "filtered_snps_" + chr + ".vcf"
            filtered_indel = use_dir + "filtered_indel_" + chr + ".vcf"
            hard_filtering = "hard_filtering_" + chr + ".vcf"
        
        
            Select_SNP_Variants_string = """
                
                
                
            echo '\\n---------- Select SNP Variants -------------\\n'
            %(GATK_path)s \\
                -T SelectVariants \\
                -R %(genome_reference)s \\
                -V %(input_file)s \\
                -selectType SNP \\
                -L %(my_chrom)s \\
                -o %(raw_snps)s
            """ % { 
                    "GATK_path" : self.params["script_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "input_file" : input_file,
                    "my_chrom" : chr,
                    "raw_snps" : raw_snps

            }      
            SNP_VariantFiltration_string = """
                
                
                
            echo '\\n---------- Create SNP VariantFiltration file -------------\\n'
            %(GATK_path)s \\
                -T VariantFiltration \\
                -R %(genome_reference)s \\
                -V %(raw_snps)s \\
                -L %(my_chrom)s \\
                --filterExpression %(filterExpression_SNP)s \\
                --filterName "my_snp_filter" \\
                -o %(filtered_snps)s
            """ % { 
                                "GATK_path" : self.params["script_path"],
                                "genome_reference" : self.params["genome_reference"],
                                "raw_snps" : raw_snps,
                                "filtered_snps" : filtered_snps,
                                "my_chrom" : chr,
                                "filterExpression_SNP" : self.params["filterExpression_SNP"]
                        }             

######################################################################### INDEL
            Select_INDEL_Variants_string = """
                
                
                
            echo '\\n---------- Select indel Variants -------------\\n'
            %(GATK_path)s \\
                -T SelectVariants \\
                -R %(genome_reference)s \\
                -V %(input_file)s \\
                -L %(my_chrom)s \\
                -selectType INDEL  \\
                -o %(raw_indel)s
            """ % { 
                                "GATK_path" : self.params["script_path"],
                                "genome_reference" : self.params["genome_reference"],
                                "input_file" : input_file,
                                "my_chrom" : chr,
                                "raw_indel" : raw_indel

                        }      
            INDEL_VariantFiltration_string = """
                
                
                
            echo '\\n---------- Create indel ApplyRecalibration file -------------\\n'
            %(GATK_path)s \\
                -T VariantFiltration \\
                -R %(genome_reference)s \\
                -V %(raw_indel)s \\
                -L %(my_chrom)s \\
                --filterExpression %(filterExpression_INDEL)s \\
                --filterName "my_indel_filter" \\
                -o %(filtered_indel)s
            """ % { 
                    "GATK_path" : self.params["script_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "raw_indel" : raw_indel,
                    "my_chrom" : chr,
                    "filterExpression_INDEL" : self.params["filterExpression_INDEL"],
                    "filtered_indel" : filtered_indel
            }


            CombineVariants_string = """
                
                
                
            echo '\\n---------- CombineVariants SNP and Indel -------------\\n'
            %(GATK_path)s \\
                -T CombineVariants \\
                -R %(genome_reference)s \\
                -V %(filtered_snps)s \\
                -V %(filtered_indel)s \\
                -L %(my_chrom)s \\
                -genotypeMergeOptions UNSORTED \\
                -o %(dir)s%(hard_filtering)s
            """ % { 
                    "GATK_path" : self.params["script_path"],
                    "genome_reference" : self.params["genome_reference"],
                    "filtered_indel" : filtered_indel,
                    "filtered_snps" : filtered_snps,
                    "dir"            : use_dir,
                    "my_chrom" : chr,
                    "hard_filtering" : hard_filtering
            }
            
            self.script = ""
        
         
            self.script = Select_SNP_Variants_string + \
                        SNP_VariantFiltration_string + \
                        Select_INDEL_Variants_string + \
                        INDEL_VariantFiltration_string + \
                        CombineVariants_string
                        
            self.sample_data[chr]["vcf"]= self.base_dir + hard_filtering
            self.stamp_file(self.sample_data[chr]["vcf"])
            
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
                    

