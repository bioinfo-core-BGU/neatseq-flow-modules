# -*- coding: UTF-8 -*-
"""
``VEP``
-----------------------------------------------------------------

:Authors: Michal Gordon
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for annotation of the multi VCF file

.. attention:: The module generates a script for each chromosome.

The programs included in the module are the following:

* ``VEP`` (`Variant Effect Predictor <https://www.ensembl.org/info/docs/tools/vep/index.html>`_. )



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``self.sample_data[chr]["vcf"]``


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``self.sample_data[chr]["vcf"]`` - annotated multi-VCF per chromosome


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "chrom_list", "list of chromosome names as mentioned in the BAM file separated by ','"

.. Note:: VEP parameters can be passed via ``redirects``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    VEP1:
        module: VEP 
        base: GATK_hard_filters1
        script_path: /path/to/vep
        chrom_list: "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT" 
        redirects:
            --format: vcf
            --offline: null
            --species: homo_sapiens
            --fork: 10
            --assembly: GRCh37
            --max_af: null
            --pick: null
            --dir: /path/to/VEP/ensembl-vep-release-88.10/cache
            --check_existing: null
            --symbol: null
            --force_overwrite: null
            --vcf: null
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
McLaren, William, et al. "The ensembl variant effect predictor." Genome biology 17.1 (2016): 122.‏

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Michal Gordon"
__version__ = "1.6.0"

class Step_VEP(Step):
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

            output_file = use_dir + self.sample_data["Title"] + "_VEP_" + chr + ".vcf"
            self.script = ""
            self.script += self.get_script_const()
            self.script += """-i %(input_file)s \\
            -o %(output_file)s\n
            """ % { 
                    "input_file" : input_file,
                    "output_file" : output_file
            }      
        

            self.sample_data[chr]["vcf"]= use_dir + self.sample_data["Title"] + "_VEP_" + chr + ".vcf"
            self.stamp_file(self.sample_data[chr]["vcf"])
            
            self.local_finish(use_dir,self.base_dir)

            self.create_low_level_script()
                        

