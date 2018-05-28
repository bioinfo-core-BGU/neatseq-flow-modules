# -*- coding: UTF-8 -*-
""" 
``trinmap_statistics``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for running ``abundance_estimates_to_matrix.pl`` on genes or isoforms counts tables produced by ``align_and_estimate_abundance.pl``

See the script documentation `here <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#building-expression-matrices>`_.

This conversion makes sense at the project level - combining all sample matrices into a single, normalized, comparison table. However, for completeness, we included a sample scope option for running the script in each sample separately.

.. Note:: ``scope`` is not defined for this module. It only makes sense to run ``abundance_estimates_to_matrix`` when comparing many samples against a single assembly

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* Either ``genes.results`` or ``isoforms.results`` files in the following slots:
        
        * ``sample_data[<sample>]["genes.results"]``
        * ``sample_data[<sample>]["isoforms.results"]``

    
Output:
~~~~~~~~~~~~~

* Creates the following files in the following slots:
    
        * ``<project>.counts.matrix``                        in ``self.sample_data["counts.matrix"]``
        * ``<project>.not_cross_norm.fpkm.tmp``              in ``self.sample_data["not_cross_norm.fpkm.tmp"]``
        * ``<project>.not_cross_norm.fpkm.tmp.TMM_info.txt`` in ``self.sample_data["not_cross_norm.fpkm.tmp.TMM_info.txt"]``
        * ``<project>.TMM.fpkm.matrix``                      in ``self.sample_data["TMM.fpkm.matrix"]``

        
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "use_genes", "", "Use 'genes.results' matrix. If not passed, use 'isoforms.results'"
    "redirects: --gene_trans_map", "path or 'none'", "If path, use path as gene_trans_map for all samples. If 'none', does not produce gene level estimates. **In order to use an internal gene_trans_map, do not pass this parameter!**"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    trin_map_stats:
        module:             trinmap_statistics
        base:               trin_map1
        script_path:        /path/to/abundance_estimates_to_matrix.pl
        use_genes:       
        redirects:
            --est_method:   RSEM

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011. **Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data**. *Nature biotechnology*, 29(7), p.644.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_trinmap_statistics(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "trin_stats"
        
        if "use_genes" not in self.params:
            self.write_warning("'use_genes' not passed. Using 'isoforms.results' matrix")

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        # In new version, --gene_trans_map is compulsory! Adding
        # If not passed:
            # If one exists, use it. 
            # Otherwise, specify "none"
        # If passed:
            # If with value, use the value and set project "gene_trans_map" to value
            # Otherwise, use existing
        if "--gene_trans_map" not in self.params["redir_params"]:
            if "gene_trans_map" in self.sample_data:
                self.params["redir_params"]["--gene_trans_map"] = self.sample_data["gene_trans_map"]
                self.use_gene_trans_map = True
            else:
                self.params["redir_params"]["--gene_trans_map"] = "none"
                self.use_gene_trans_map = False
                    
        else:  # --gene_trans_map is defined in redir_params
            if self.params["redir_params"]["--gene_trans_map"] == None:
                raise AssertionExcept("You passed --gene_trans_map with no value. Please specify path or 'none'")
            elif self.params["redir_params"]["--gene_trans_map"] == "none":
                self.use_gene_trans_map = False
            else:
                self.sample_data["gene_trans_map"] = self.params["redir_params"]["--gene_trans_map"]
                self.use_gene_trans_map = True
            
                
            
        
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        

    def build_scripts(self):
        

        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        
        prefix = self.sample_data["Title"]

        self.script += self.get_script_const()

            
        self.script += "--out_prefix %s \\\n\t" % os.sep.join([use_dir, prefix])
        # type2use is 'genes.results' or 'isoforms.results'. This is used to then select the correct slot from "mapping"
        type2use = "genes.results" if "use_genes" in self.params.keys() else "isoforms.results"
        
        for sample in self.sample_data["samples"]:
            try:
                self.script += "%s \\\n\t" % self.sample_data[sample][type2use] 
            except:
                raise AssertionExcept("file type %s does not exist for sample." % type2use, sample)
        
        self.script = self.script.rstrip("\\\n\t")
        self.script += "\n\n"
        
        
        if not "version" in self.params or self.params["version"].lower() == "new":

        # Storing all output files even though probably not very useful downstream...
            self.sample_data["isoform.raw_counts"] = os.sep.join([self.base_dir,  "%s.isoform.counts.matrix" % prefix])
            self.sample_data["isoform.norm_counts"] = os.sep.join([self.base_dir, "%s.isoform.TPM.not_cross_norm" % prefix])
            
            self.stamp_file(self.sample_data["isoform.raw_counts"] )
            self.stamp_file(self.sample_data["isoform.norm_counts"])
            
            
            if(self.use_gene_trans_map):  # True when --gene_trans_map is not "none"
                self.sample_data["gene.raw_counts"] = os.sep.join([self.base_dir,  "%s.gene.counts.matrix" % prefix])
                self.sample_data["gene.norm_counts"] = os.sep.join([self.base_dir, "%s.gene.TPM.not_cross_norm" % prefix])
                self.stamp_file(self.sample_data["gene.raw_counts"] )
                self.stamp_file(self.sample_data["gene.norm_counts"])
            
        
        else:
            self.write_warning("Not storing output files for old version of trinity. If required, load the appropriate files with a 'manage_types' module")

       
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     
            
        
        
        self.create_low_level_script()
                    
    
     