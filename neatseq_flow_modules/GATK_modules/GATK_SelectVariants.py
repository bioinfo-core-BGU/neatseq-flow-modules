# -*- coding: UTF-8 -*-
"""
``GATK_SelectVariants``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for separation of multi-VCF per-chromosome to one VCF per-sample per-chromosome

.. attention:: The module generates a script for each sample/chromosome.

The programs included in the module are the following:

* ``SelectVariants`` (GATK) 



Requires
~~~~~~~~~~~~


* ``self.sample_data[chr]["vcf"]``
* ``self.params["genome_reference"]``
* ``self.params["chrom_list"]`` - list of chromosomes names as mentioned in BAM file separated by ','


Output
~~~~~~~~~~~~~

* ``self.sample_data[sample][chr]["GATK_vcf"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "genome_reference", "path to reference genome", "..."
    "chrom_list", "", "list of chromosome names as mentioned in the BAM file separated by ',' "

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    GATK_SelectVariants_VEPfiltered:
        module: GATK_SelectVariants
        base: VEP1
        script_path: /path/to/GenomeAnalysisTK.jar        
        chrom_list: "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT" 
        genome_reference:   /path/to/gatk/bundle/b37/human_g1k_v37_decoy.fasta
        redirects:
            --setFilteredGtToNocall: null

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Van der Auwera, Geraldine A., et al. "From FastQ data to high‐confidence variant calls: the genome analysis toolkit best practices pipeline." Current protocols in bioinformatics 43.1 (2013): 11-10.‏

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

class Step_GATK_SelectVariants (Step):
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
        for chr in self.params["chrom_list"].split(','):
            chr = chr.strip()

            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
                # Name of specific script:
                self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample,chr])
                self.script = ""
                
                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
                self.script += self.get_script_const()
                
                output_file = sample_dir + sample + "_" + chr + "_GATK.vcf" 
                
                my_string = """-T SelectVariants \\
        -R %(genome_reference)s \\
        -V %(input_full_vcf)s \\
        -L %(my_chrom)s \\
        -o %(output_sample_vcf)s \\
        -sn %(SAMPLE)s

                """ % { "sample_dir" : sample_dir,
                        "genome_reference" : self.params["genome_reference"],
                        "input_full_vcf" : self.sample_data[chr]["vcf"],
                        "my_chrom" : chr,
                        "output_sample_vcf" : output_file,
                        "SAMPLE" : sample
                }
                final_output = sample_dir + sample + "_" + chr + "_GATK_final.vcf" 
# delete as reference rows                 
                my_string_egrep = """egrep -v "AC=0;|AC=0,0;|AC=0,0,0;|AC=0,0,0,0;|AC=0,0,0,0,0;|AC=0,0,0,0,0,0;"  """ + output_file + " > " + final_output
                my_rm_script = "rm -f " +  output_file + "\n\t" + "rm -f " +  output_file + ".idx"
                self.script += my_string + "\n\t" + my_string_egrep + "\n\t" + my_rm_script
                #self.get_script_env_path()
                

                self.sample_data[sample][chr] = {}
                self.sample_data[sample][chr]["GATK_vcf"] = final_output
                self.stamp_file(self.sample_data[sample][chr]["GATK_vcf"])
                
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

                self.create_low_level_script()
                        
