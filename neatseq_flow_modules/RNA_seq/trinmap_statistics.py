# -*- coding: UTF-8 -*-
""" 
``trinity``
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

    "use_isoforms", "", "Use 'isoforms.results' matrix. If not passed, use 'genes.results'"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    trin_map_stats:
        module:             trinmap_statistics
        base:               trin_map1
        script_path:        /path/to/abundance_estimates_to_matrix.pl
        use_isoforms:       
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
        
        if "use_isoforms" not in self.params:
            self.write_warning("'use_isoforms' not passed. Using 'genes.results' matrix")
        # if self.params["use_isoforms"] not in ['genes.results' , 'isoforms.results']:
            # raise AssertionExcept("'use_isoforms' can be either 'genes.results' or 'isoforms.results'")
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        
        
        
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
        # In new version, --gene_trans_map is compulsory! Adding
        if "--gene_trans_map" not in self.params["redir_params"]:
            if "gene_trans_map" not in self.sample_data:
                self.script += "--gene_trans_map none \\\n\t"
            else:
                self.script += "--gene_trans_map %s \\\n\t" % self.sample_data["gene_trans_map"]

        self.script += "--out_prefix %s \\\n\t" % os.sep.join([use_dir, prefix])
        # type2use is 'genes.results' or 'isoforms.results'. This is used to then select the correct slot from "mapping"
        type2use = "isoforms.results" if "use_isoforms" in self.params.keys() else "genes.results"
        
        for sample in self.sample_data["samples"]:
            try:
                self.script += "%s \\\n\t" % self.sample_data[sample][type2use] 
            except:
                raise AssertionExcept("file type %s does not exist for sample." % type2use, sample)
        
        self.script.rstrip("\\\n\t")
        self.script += "\n\n"
        
        
        # Storing all output files even though probably not very useful downstream...
        self.sample_data["counts.matrix"] = os.sep.join([self.base_dir, "%s.counts.matrix" % prefix])
        self.sample_data["not_cross_norm.fpkm.tmp"] = os.sep.join([self.base_dir, "%s.not_cross_norm.fpkm.tmp" % prefix])
        self.sample_data["not_cross_norm.fpkm.tmp.TMM_info.txt"] = os.sep.join([self.base_dir, "%s.not_cross_norm.fpkm.tmp.TMM_info.txt" % prefix])
        self.sample_data["TMM.fpkm.matrix"] = os.sep.join([self.base_dir, "%s.TMM.fpkm.matrix" % prefix])

        self.stamp_file(self.sample_data["counts.matrix"])
        self.stamp_file(self.sample_data["not_cross_norm.fpkm.tmp"])
        self.stamp_file(self.sample_data["not_cross_norm.fpkm.tmp.TMM_info.txt"])
        self.stamp_file(self.sample_data["TMM.fpkm.matrix"])


       
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     
            
        
        
        self.create_low_level_script()
                    
    
     