# -*- coding: UTF-8 -*-
""" 
``Trinity_gene_to_trans_map`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for RNA_seq assembly using the `Trinity assembler`_.

.. _Trinity assembler: https://github.com/trinityrnaseq/trinityrnaseq/wiki
 
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    * ``fastq`` files in at least one of the following slots:
        
        * ``sample_data[<sample>]["fastq.F"]``
        * ``sample_data[<sample>]["fastq.R"]``
        * ``sample_data[<sample>]["fastq.S"]``

    
Output:
~~~~~~~~~~~~~

    * puts ``fasta`` output files in the following slots:
        
        * for sample-wise assembly:
        
            * ``sample_data[<sample>]["fasta.nucl"]``
            * ``sample_data[<sample>]["Trinity.contigs"]``
        
        * for project-wise assembly:
        
            * ``sample_data["fasta.nucl"]``
            * ``sample_data["Trinity.contigs"]``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Create one assembly for all samples or one assembly per sample."
    "skip_gene_to_trans_map", "", "Set to skip executing get_Trinity_gene_to_trans_map.pl on assembly"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trinity1:
        module:     trinity
        base:       trin_tags1
        script_path: /path/to/Trinity
        qsub_params:
            node:      sge213
            -pe:       shared 20
        # skip_gene_to_trans_map:
        redirects:
            --grid_conf:        /path/to/SGE_Trinity_conf.txt
            --CPU:              20
            --seqType:          fq
            --JM:               140G
            --min_kmer_cov:     2
            --full_cleanup:

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011. **Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data**. *Nature biotechnology*, 29(7), p.644.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_Trinity_gene_to_trans_map(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        if "trinity" not in [self.pipe_data["names_index"][step] for step in self.get_depend_list()]:
            self.write_warning("No trinity in history. Are you sure of what you are attempting to do?")
        
        
        
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                if not "fasta.nucl" in self.sample_data:
                    raise AssertionExcept("No fasta file of type 'nucl' in project\n")

            elif self.params["scope"]=="sample":
                
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    if not "fasta.nucl" in self.sample_data[sample]:
                        raise AssertionExcept("No fasta file of type 'nucl'\n",sample)
                    
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
        
        
        ##########################
        

            
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):
    
        if self.params["scope"] == "project":
            self.build_scripts_project()
        else:
            self.build_scripts_sample()
            
            
    def build_scripts_project(self):
        
        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        output_basefn = "%s.gene_trans_map" % os.path.basename(self.sample_data["fasta.nucl"])
        
        self.script = """
{script_path}{transcriptome} \\
	> {map} 
""".format(transcriptome = "%s" % self.sample_data["fasta.nucl"],
           script_path   = self.get_script_const(),
           map           = os.path.join(use_dir, output_basefn))
           
        
        
        self.sample_data["gene_trans_map"] = os.path.join(self.base_dir, output_basefn)
        self.stamp_file(self.sample_data["gene_trans_map"])
       
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        
        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            output_basefn = "%s.gene_trans_map" % os.path.basename(self.sample_data[sample]["fasta.nucl"])
            
            self.script = """
{script_path}{transcriptome} \\
	> {map} 
""".format(transcriptome = "%s" % self.sample_data[sample]["fasta.nucl"],
           script_path   = self.get_script_const(),
           map           = os.path.join(use_dir, output_basefn))
               
            
            
            self.sample_data[sample]["gene_trans_map"] = os.path.join(sample_dir, output_basefn)
            self.stamp_file(self.sample_data[sample]["gene_trans_map"])
           
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                        
            
            
                 
            
     