# -*- coding: UTF-8 -*-
""" 
``merge_project``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for merging sample `fasta` files into a single project wide `fasta` file.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A `fasta` file in one of the following slots (scope = "sample"):

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

* A `fasta` file in one of the following slots (scope = "project"):

    * ``sample_data["fasta.nucl"]``
    * ``sample_data["fasta.prot"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slots (scope = "sample"):
        
    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

* Puts output files in the following slots (scope = "project"):

    * ``sample_data["fasta.nucl"]``
    * ``sample_data["fasta.prot"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    split_fasta1:
        module:         split_fasta
        base:           merge1
        script_path:    cat
        type:           fasta.nucl


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept
from astropy.vo import samp


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_split_fasta(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        if "type" not in self.params:
            raise AssertionExcept("Please supply a 'type' parameter: 'nucl' or 'prot'")
        if self.params["type"] not in ["nucl","prot"]:
                raise AssertionExcept("'type' parameter must be 'nucl' or 'prot'")
        if "scope" not in self.params:
            raise AssertionExcept("Please supply a 'scope' parameter: 'sample' or 'project'")
        if self.params["scope"] not in ["sample","project"]:
                raise AssertionExcept("'scope' parameter must be 'sample' or 'project'")

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Check that scope exists and fasta.nucl or fasta.prot exist
        # Check that subsample_num exists
        # Check that type exists

        self.params["type"] = "fasta.{type}".format(type=self.params["type"])
        if self.params["scope"]=="project":
            if self.params["type"] not in self.sample_data:
                raise AssertionExcept("{type} does not exist in project".format(type=self.params["type"]))
        else:
            for sample in self.sample_data["samples"]:
                if self.params["type"] not in self.sample_data[sample]:
                    raise AssertionExcept("{type} does not exist in sample".format(type=self.params["type"]), sample)
                
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        if self.params["scope"] == "project":
            
                            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            
            self.script = """
FASTA={project_fasta}
SEQPERFRAG=$[$(grep -c "^>" $FASTA)/$[{subsample_num}-1]]
awk -v seqs="$SEQPERFRAG" 'BEGIN {{n_seq=0; file_cnt=1;}} /^>/ {{ if(n_seq%seqs==0){{file=sprintf("{use_dir}subsample%d.fa",file_cnt); file_cnt++; }} print >> file; n_seq++; next;}} {{ print >> file; }}' < $FASTA

""".format(project_fasta = self.sample_data[self.params["type"]],
            subsample_num = self.params["subsample_num"],
            use_dir = use_dir)

            self.sample_data["samples.old"] = self.sample_data["samples"]            
            self.sample_data["samples"] = ["subsample{num}".format(num=num) for num in range(1,self.params["subsample_num"]+1)]
                
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.sample_data[sample] = dict()
                self.sample_data[sample][self.params["type"]] = "{use_dir}{sample}.fa".format(use_dir=self.base_dir,sample=sample)
                self.stamp_file(self.sample_data[sample][self.params["type"]])
                

                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                    
        
        else: # self.params["scope"] == "sample"
        
            raise AssertionExcept("Not defined yet...")
            # Each iteration must define the following class variables:
                # spec_script_name
                # script
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                    
                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                
                self.create_low_level_script()
                    