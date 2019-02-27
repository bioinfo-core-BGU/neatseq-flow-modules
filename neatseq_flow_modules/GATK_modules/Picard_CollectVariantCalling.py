# -*- coding: UTF-8 -*-
"""
``Picard_CollectVariantCalling``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for generating SNP and indel statistics information


The programs included in the module are the following:

* ``CollectVariantCallingMetrics`` Picard tool to generate A collection of metrics relating to snps and indels within a variant-calling file (VCF) 


Requires
~~~~~~~~~~~~

* A fastq file in the following location:

    * ``self.sample_data[chr]["vcf"]``

Output
~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Picard_CollectVariantCalling1:
        module: Picard_CollectVariantCalling 
        base: GATK_hard_filters1
        script_path: /path/to/java -jar /path/to/picard.jar
        DBSNP: /path/to/bundle/b37/dbsnp_138.b37.vcf
        chrom_list: "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT"




References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
http://broadinstitute.github.io/picard/

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"

class Step_Picard_CollectVariantCalling(Step):
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
            self.script = ""
            chr = chr.strip()
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],chr])

            input_file = self.sample_data[chr]["vcf"]

        
            DBSNP = self.params["DBSNP"]
            output_CollectVariantCallingMetrics =  use_dir + self.sample_data["Title"] + "CollectVariantCallingMetrics_" + chr + ".txt" 
            my_string = """
                    cd %(dir)s
                    echo '\\n---------- CollectVariantCallingMetrics -------------\\n'
                    %(picard_path)s CollectVariantCallingMetrics \\
                    INPUT=%(input_file)s \\
                    DBSNP=%(DBSNP)s \\
                    OUTPUT=%(output_file)s

                """ % { 
                        "picard_path" : self.params["script_path"],
                        "input_file" : input_file,
                        "DBSNP" : DBSNP,
                        "dir" : use_dir,
                        "output_file" : output_CollectVariantCallingMetrics
                }
                
                
                
            self.script += my_string
            #self.get_script_env_path()
                
                                       
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()

