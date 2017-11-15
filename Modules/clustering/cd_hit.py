# -*- coding: UTF-8 -*-
""" 
``cd_hit``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for clustering with cd-hit/ch-hit-est:

This module runs both cd-hit and cd-hit-est. The type of sequence (`nucl` or `prot`) will be determined by the program supplied in `script_path`.

You must make sure that the required file exists: If clustering prot sequences with ``cd-hit-est``, make sure there is a ``fasta.prot`` file, etc.

**CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences**, Weizhong Li & Adam Godzik. *Bioinformatics*, (2006) 22:1658-1659

**CD-HIT: accelerated for clustering the next generation sequencing data**, Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. *Bioinformatics*, (2012) 28:3150-3152


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fasta files in the following slot (scope = sample):

    * ``sample_data[<sample>]["fasta.nucl"|"fasta.prot"]``
    
* fasta files in the following slot (scope = project):

    * ``sample_data["fasta.nucl"|"fasta.prot"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts the output fasta file in the fasta slot:

    ``self.sample_data[<sample>]["fasta.nucl"|"fasta.prot"]``

* Or

    ``self.sample_data["fasta.nucl"|"fasta.prot"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Indicates whether to use a project or sample fasta."


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::

    clust_proj:
        module: cd_hit
        base: derepel_proj
        script_path: 'path/to/cd-hit-est'
        qsub_params:
            -pe: shared 40
        scope: project
        redirects:
            -T: 40

            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fu, L., Niu, B., Zhu, Z., Wu, S. and Li, W., 2012. **CD-HIT: accelerated for clustering the next-generation sequencing data.** *Bioinformatics*, 28(23), pp.3150-3152.
"""


import os, sys, re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_cd_hit(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
    
        if "scope" not in self.params:
            raise AssertionExcept("You must specify 'scope'\n")

        if re.search("cd-hit-est$", self.params["script_path"]):
            self.type = "nucl"
        elif re.search("cd-hit$", self.params["script_path"]):
            self.type = "prot"
        else:
            raise AssertionExcept("'script_path' must include either 'cd-hit' or 'cd-hit-est'. You set it to '%s'" % self.params["script_path"])
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if self.params["scope"] == "sample":
            pass
            # raise AssertionExcept("sample scope not yet defined... Sorry")
        elif self.params["scope"] == "project":
            pass

        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

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
        
        
        # Each iteration must define the following class variables:
            # self.spec_script_name
            # self.script
        if self.params["scope"] == "project":

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            # output_prefix = sample + "_bowtie2_map"

            if self.type == "nucl":
                try:
                    input_file = self.sample_data["fasta.nucl"]
                except:
                    raise AssertionExcept("`nucl` fasta file does not exist at project scope. Did you mean cd-hit instead of cd-hit-est?")
            else: # == "prot"
                try:
                    input_file = self.sample_data["fasta.prot"]
                except:
                    raise AssertionExcept("`prot` fasta file does not exist at project scope. Did you mean cd-hit-est instead of cd-hit?")
            
            output_prefix = os.path.basename(input_file)
            
            # Get constant part of script:
            self.script += self.get_script_const()
            self.script += "-i {infn} \\\n\t".format(infn = input_file)
            self.script += "-o {outdir}{ossep}{outfn} \n\n".format(outdir=use_dir, ossep = os.sep, outfn = output_prefix)

            
            self.sample_data["fasta." + self.type] = "{outdir}{ossep}{outfn}".format(outdir = self.base_dir, ossep = os.sep, outfn = output_prefix)
            self.sample_data["cd_hit." + self.type] = self.sample_data["fasta." + self.type]
            self.stamp_file(self.sample_data["fasta." + self.type])
                    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                            
        else:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
     
                # Define location and prefix for output files:
                # output_prefix = sample + "_bowtie2_map"


                if self.type == "nucl":
                    try:
                        input_file = self.sample_data[sample]["fasta.nucl"]
                    except:
                        raise AssertionExcept("`nucl` fasta file does not exist at project scope. Did you mean cd-hit instead of cd-hit-est?")
                else: # == "prot"
                    try:
                        input_file = self.sample_data[sample]["fasta.prot"]
                    except:
                        raise AssertionExcept("`prot` fasta file does not exist at project scope. Did you mean cd-hit-est instead of cd-hit?")
                                
                output_prefix = os.path.basename(input_file)
                
                
                # Get constant part of script:
                self.script += self.get_script_const()
                self.script += "-i {infn} \\\n\t".format(infn = input_file)
                self.script += "-o {outdir}{ossep}{outfn} \n\n".format(outdir=use_dir, ossep = os.sep, outfn = output_prefix)

                
                self.sample_data[sample]["fasta." + self.type] = "{outdir}{ossep}{outfn}".format(outdir = sample_dir, ossep = os.sep, outfn = output_prefix)
                self.sample_data[sample]["cd_hit." + self.type] = self.sample_data[sample]["fasta." + self.type]
                self.stamp_file(self.sample_data[sample]["fasta." + self.type])
                        
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                              
                
