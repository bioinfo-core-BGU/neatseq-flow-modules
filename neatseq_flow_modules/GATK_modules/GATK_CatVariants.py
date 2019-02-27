# -*- coding: UTF-8 -*-
"""
``GATK_CatVariants``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module to concatenate chromosome to get one VCF file for each sample.

.. attention:: The module generate script for each sample - chromosom.

The programs included in the module are the following:

* ``CatVariants`` (GATK) 



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``self.sample_data[sample][chr]["GATK_vcf"]``


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``self.sample_data[sample]["vcf"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "genome_reference", "", ""
    "chrom_list", "Comma-separated list of chromosome names as mentioned in the BAM file"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    GATK_CatVariants1:
        module: GATK_CatVariants
        base: GATK_SelectVariants_VEPfiltered
        script_path:     /path/to/java -cp /path/to/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants
        genome_reference:   /path/to/gatk/bundle/b37/human_g1k_v37_decoy.fasta
        chrom_list: "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT"


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Van der Auwera, Geraldine A., et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43.1 (2013): 11-10.‏

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"


class Step_GATK_CatVariants(Step):
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
            output_file = self.base_dir + sample + "_CatVariants.vcf"
            my_CatVariants_string = """
                    
                    
                    
    echo '\\n---------- Create SNP VariantRecalibrator file -------------\\n'
    %(GATK_path)s \\
        -R %(genome_reference)s \\
    """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"]
                }      
            for chr in self.params["chrom_list"].split(','):
                chr = chr.strip()

                # Name of specific script:
                my_CatVariants_string = my_CatVariants_string + "    -V " + self.sample_data[sample][chr]["GATK_vcf"] + " \\\n"
            my_CatVariants_string = my_CatVariants_string + "    -out " + output_file + " \n"
            self.script = my_CatVariants_string
    #        print self.script
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample])
            # self.spec_script_name = set_spec_script_name()
            # self.jid_name_sep instead of "_"
            use_dir = self.local_start(self.base_dir)
            self.sample_data[sample]["vcf"] = output_file
            self.stamp_file(self.sample_data[sample]["vcf"])
                
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
    '''
     java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
        -R reference.fasta \
        -V input1.vcf \
        -V input2.vcf \
        -out output.vcf \
        -assumeSorted'''

