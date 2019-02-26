# -*- coding: UTF-8 -*-
"""
``GATK_merge_gvcf``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for combine g.vcf files to cohorts.

.. attention:: The module generate script for each sample-chromosom.

The programs included in the module are the following:

* ``CombineGVCFs`` (GATK) 



**Requires**:


    * ``self.params["cohort_size"]`` - number of g.vcf file to be in each cohort
    * ``self.params["genome_reference"]``
    * ``self.params["chrom_list"]`` - list of chromosomes names as mentioned in BAM file separated by ','
    * ``self.sample_data[sample][chr]["GATK_g.vcf"]``

**Output**:

    * ``self.sample_data["cohorts"]``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    gatk_merge_gvcf:
        module: GATK_merge_gvcf
        base: GATK_gvcf
        script_path:     /path/to/java -jar /path/to/GenomeAnalysisTK.jar
        genome_reference:    /path/to/gatk/bundle/b37/human_g1k_v37_decoy.fasta
        cohort_size: 10
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

class Step_GATK_merge_gvcf(Step):
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
            
    def chunks(self, l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]      
    
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
    
        my_cohorts = self.chunks(self.sample_data["samples"],self.params["cohort_size"])
        my_cohort_counter = 0
        

        new_sample_list = list()
        # print self.params["chrom_list"]
        for cohort in my_cohorts:
            # print cohort

            my_cohort_counter = my_cohort_counter +1
            cohort_name = "cohort{n}".format(n = my_cohort_counter)
            self.sample_data[cohort_name] = {}
                
            new_sample_list.append(cohort_name)
            
            for chr in self.params["chrom_list"].split(','):
                chr = chr.strip()
                my_variant_string = ""

                
                for sample in cohort:
                    my_variant_string += "\t--variant " + self.sample_data[sample][chr]["GATK_g.vcf"] + " \\\n"     # Getting list of samples within cohort
#                   self.sample_data[sample][chr]["cohort"] = cohort_name
                
                
                # Name of specific script:
                self.spec_script_name = self.jid_name_sep.join([self.step,self.name,cohort_name,"chr",chr])
                self.script = ""
                
                # Make a dir for the current cohort:
                sample_dir = self.make_folder_for_sample(cohort_name)
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
                
                #cd %(sample_dir)s
                my_string = """
echo '\\n---------- Create gvcf file -------------\\n'
%(GATK_path)s \\
    -T CombineGVCFs \\
    -R %(genome_reference)s \\
%(my_variant)s\t-o %(output_cohort_gvcf_creation)s
                """ % { 
                        "GATK_path" : self.params["script_path"],
                        "genome_reference" : self.params["genome_reference"],
                        "output_cohort_gvcf_creation" : "{d}{n}.chr{c}.g.vcf".format(d = sample_dir, n = cohort_name, c= chr),
                        "my_variant" : my_variant_string                    
                }
                self.script += my_string
                #self.get_script_env_path()
                
                self.sample_data[cohort_name][chr] = {}
                self.sample_data[cohort_name][chr]["g.vcf"] = "{d}{n}.chr{c}.g.vcf".format(d = sample_dir, n = cohort_name, c= chr)
                self.stamp_file(self.sample_data[cohort_name][chr]["g.vcf"])
                
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

                self.create_low_level_script()
                        
            self.sample_data["cohorts"] = new_sample_list  #???????
        
#        print "From merge_gvcf::\n"
#        print self.sample_data["cohorts"][chr]
