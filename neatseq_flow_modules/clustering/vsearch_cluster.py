# -*- coding: UTF-8 -*-
""" 
``vsearch_cluster``
--------------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running vsearch clustering:

The reads stored in fasta files are clustered with one of the 3 methods available: cluster_fast, cluster_size or cluster_smallmem.

..Note: At the moment this works on the `nucl` fasta only. See the web: https://github.com/torognes/vsearch/issues/42

Output types are defined with the `outputs` parameter which can be a comma separated list of the following:

    biomout,mothur_shared_out,otutabout,profile,uc

Fasta output files are defined with the `fasta_outputs` parameter which can be a comma separated list of the following:

    centroids,consout,msaout

By default, the `centroids` file is stored in the `fasta` slot. Change this by setting `store_fasta` to one of the types listed above, i.e. centroids,consout or msaout

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fasta files in the following slot (scope = sample):

    * ``sample_data[<sample>]["fasta.nucl"]``
    
* fasta files in the following slot (scope = project):

    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts required output in similarly named slots, e.g.:

    ``self.sample_data[<sample>]["vsearch.centroids"]`` or
    ``self.sample_data["project_data"]["vsearch.centroids"]``


* Puts the required fasta in the fasta slot:

    ``self.sample_data[<sample>]["fasta.nucl"]`` or ``self.sample_data["project_data"]["fasta.nucl"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "outputs", "biomout,mothur_shared_out,otutabout,profile,uc", "List of outputs other than fasta type outputs (see `fasta_outputs`"
    "fasta_outputs", "centroids,consout,msaout", "A list of fasta types to produce."
    "store_fasta", "centroids|consout|msaout", "The fasta type to store in fasta slot"
    "scope", "project | sample", "Indicates whether to use a project or sample nucl fasta."



Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::

    clust_proj:
        module: vsearch_cluster
        base: derepel_proj
        script_path: '{Vars.vsearch_path}/vsearch'
        qsub_params:
            -pe: shared 40
        fasta_outputs: centroids,consout
        outputs: uc
        store_fasta: centroids
        scope: project
        type: cluster_fast
        redirects:
            --id: 0.85  # From ipyrad defaults
            --qmask: dust
            --strand: both
            --threads: 40
            --sizein:
            --sizeout:

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rognes, T., Flouri, T., Nichols, B., Quince, C. and Mahé, F., 2016. **VSEARCH: a versatile open source tool for metagenomics**. *PeerJ*, 4, p.e2584.

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_vsearch_cluster(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
    
        if "scope" not in self.params:
            raise AssertionExcept("You must specify 'scope'\n")
        if "type" not in self.params:
            raise AssertionExcept("You must specify 'type': Either cluster_fast or cluster_size or cluster_smallmem")
        if self.params["type"] not in ["cluster_fast","cluster_size","cluster_smallmem"]:
            raise AssertionExcept("You must specify 'type': Either cluster_fast or cluster_size or cluster_smallmem")

        # Splitting outputs list
        if "outputs" in list(self.params.keys()):
            self.params["outputs"] = self.params["outputs"].split(",")
        # Splitting fasta_outputs list
        if "fasta_outputs" in list(self.params.keys()):
            self.params["fasta_outputs"] = self.params["fasta_outputs"].split(",")
        else: # If the list isn't passed, use centroids by default
            self.params["fasta_outputs"] = ["centroids"]
            self.write_warning("Outputting centroids fasta")
        if "store_fasta" in list(self.params.keys()):
            if self.params["store_fasta"] not in self.params["fasta_outputs"]:
                raise AssertionExcept("'store_fasta' must also appear in the 'fasta_outputs' list.")
        else:
            self.params["store_fasta"] = ["centroids"]
            self.write_warning("Using centroids as default fasta")
            
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if self.params["scope"] == "sample":
        # Initializing a "fasta" dict for each sample:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have fasta nucl data" , sample)
        elif self.params["scope"] == "project":
            try:
                self.sample_data["project_data"]["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("No project-wide fasta nucl data" )

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
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
 
            # Define location and prefix for output files:
            # output_prefix = sample + "_bowtie2_map"

            input_file = self.sample_data["project_data"]["fasta.nucl"]
            
            output_prefix = os.path.basename(input_file)
            
            # Get constant part of script:
            self.script += self.get_script_const()
            if "outputs" in self.params:
                for output in self.params["outputs"]:
                    self.script += "--%s %s_%s \\\n\t" % (output, use_dir + output_prefix, output)
            for output in self.params["fasta_outputs"]:
                self.script += "--%s %s_%s.fasta \\\n\t" % (output, use_dir + output_prefix, output)
            if "centroids" not in self.params["fasta_outputs"]:
                self.script += "--centroids %s_centroids.fasta \\\n\t" % (use_dir + output_prefix)
            self.script += "--%s %s \n\n" % (self.params["type"], input_file)
            

            self.sample_data["project_data"]["fasta.nucl"] = "%s_centroids.fasta" % (self.base_dir + output_prefix)
            self.sample_data["project_data"]["vsearch.centroids"] = self.sample_data["project_data"]["fasta.nucl"]
            self.stamp_file(self.sample_data["project_data"]["fasta.nucl"])
                    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
                            
        else:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
     
                # Define location and prefix for output files:
                # output_prefix = sample + "_bowtie2_map"

                input_file = self.sample_data[sample]["fasta.nucl"]
                
                output_prefix = os.path.basename(input_file)
                
                # Get constant part of script:
                self.script += self.get_script_const()
                if "outputs" in self.params:
                    for output in self.params["outputs"].split(","):
                        self.script += "--%s %s_%s \\\n\t" % (output, use_dir + output_prefix, output)
                for output in self.params["fasta_outputs"].split(","):
                    self.script += "--%s %s_%s.fasta \\\n\t" % (output, use_dir + output_prefix, output)
                if "centroids" not in self.params["fasta_outputs"].split(","):
                    self.script += "--centroids %s_centroids.fasta \\\n\t" % (use_dir + output_prefix)
                self.script += "--%s %s \n\n" % (self.params["type"], input_file)
                

                self.sample_data[sample]["fasta.nucl"] = "%s_centroids.fasta" % (sample_dir + output_prefix)
                self.sample_data[sample]["vsearch.centroids"] = self.sample_data[sample]["fasta.nucl"]
                self.stamp_file(self.sample_data[sample]["fasta.nucl"])
        
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
           
                
                
                self.create_low_level_script()
                    
