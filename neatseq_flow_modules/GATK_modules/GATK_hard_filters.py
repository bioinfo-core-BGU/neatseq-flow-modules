# -*- coding: UTF-8 -*-
"""
``GATK_hard_filters``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for apply hard filters to a variant callset that is too small for VQSR or for which truth/training sets are not available..

.. attention:: The module generate script for each chromosom.

The programs included in the module are the following:

* ``SelectVariants and VariantFiltration`` (GATK) 



**Requires**:


    * ``self.sample_data[chr]["vcf"]``
    * ``self.params["genome_reference"]``
    * ``self.params["chrom_list"]`` - list of chromosomes names as mentioned in BAM file separated by ','
    * ``self.params["filterExpression_SNP"]`` - filter e xpression for SNP
    * ``self.params["filterExpression_INDEL"]`` - filter e xpression for INDEL


**Output**:

    * ``self.sample_data[chr]["vcf"]``


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    GATK_hard_filters1:
        module: GATK_hard_filters 
        base: GenotypeGVCFs1
        script_path:     /path/to/java -jar /path/to/GenomeAnalysisTK.jar
        genome_reference:   /path/to/gatk/bundle/b37/human_g1k_v37_decoy.fasta
        chrom_list: "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT" 
        filterExpression_SNP: '"QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'
        filterExpression_INDEL: '"QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0 || InbreedingCoeff < -0.8"'


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Van der Auwera, Geraldine A., et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43.1 (2013): 11-10.‏

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

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
                    

