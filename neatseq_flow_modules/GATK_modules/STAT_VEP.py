# -*- coding: UTF-8 -*-

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

class Step_STAT_VEP(Step):
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

        
       
        # Each iteration must define the following class variables:
        # spec_qsub_name
        # spec_script_name
        # script
######################################################## SNP

        for sample in self.sample_data["samples"]:
            output_file = self.base_dir + sample + "_stat_vep.txt"
            my_statVEP_string = """
                 
    echo '\\n---------- Create SNP statistics -------------\\n'
 
    """       
            my_printf_stop = 'printf "stop_gain: " > ' + output_file + "\n"
            my_stop_gain = ' egrep "PASS.*stop_gained" ' + self.sample_data[sample]["vcf"] + " | wc -l >> " + output_file + "\n"
            my_printf_framshift = 'printf "framshift: " >> ' + output_file + "\n"
            my_framshift = ' egrep "PASS.*frameshift" ' + self.sample_data[sample]["vcf"] + " | wc -l >> " + output_file + "\n"
            my_printf_missense_varian = 'printf "my_missense_varian: " >> ' + output_file + "\n"
            my_missense_varian = ' egrep "PASS.*missense_varian" ' + self.sample_data[sample]["vcf"] + " | wc -l >> " + output_file + "\n"
            my_printf_synonymous_variant = 'printf "my_synonymous_variant: " >> ' + output_file + "\n"
            my_synonymous_variant = ' egrep "PASS.*synonymous_variant" ' + self.sample_data[sample]["vcf"] + " | wc -l >> " + output_file + "\n"

            
            self.script = my_statVEP_string  + my_printf_stop + my_stop_gain + my_printf_framshift + my_framshift + my_printf_missense_varian +  my_missense_varian + my_printf_synonymous_variant + my_synonymous_variant
    #        print self.script
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample])
            # self.spec_script_name = set_spec_script_name()
            # self.jid_name_sep instead of "_"
            use_dir = self.local_start(self.base_dir)
            self.sample_data[sample]["vcf_stat"] = output_file
            self.stamp_file(self.sample_data[sample]["vcf_stat"])
                
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
    '''
     java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
        -R reference.fasta \
        -V input1.vcf \
        -V input2.vcf \
        -out output.vcf \
        -assumeSorted'''

