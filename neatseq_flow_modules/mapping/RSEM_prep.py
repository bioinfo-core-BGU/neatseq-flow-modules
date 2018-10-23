# -*- coding: UTF-8 -*-
""" 
``RSEM_prep``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``rsem-prepare-reference``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``  (``scope`` = ``project``)
    * ``sample_data[<sample>]["fasta.nucl"]``   (``scope`` = ``sample``)
    
* If neither exists, please supply ``reference`` parameter.
       

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Puts output index files in one of the following slot:

    * ``self.sample_data[<sample>]["RSEM_index"]``
    * ``self.sample_data["project_data"]["RSEM_index"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Where to take the reference from"
    "reference", "path to reference", "Use this fasta file. See the definition for reference_fasta_file(s) in the ARGUMENTS section of rsem-prepare-reference help"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::


    RSEM_prep_ind:
        module:             RSEM_prep
        base:               merge1
        script_path:        /path/to/RSEM
        reference:              /path/to/fasta
        redir_params:
            --gtf:          /path/to/gtf
            --transcript-to-gene-map: /path/to/map_file
    
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_RSEM_prep(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        if "reference" not in self.params:
            if "scope" not in self.params:
                raise AssertionExcept("Please supply a scope parameter: either 'sample' or 'project'!")
            elif self.params["scope"]=="sample":
                for sample in self.sample_data["samples"]:
                    if "fasta.nucl" not in self.sample_data[sample]:
                        raise AssertionExcept("No fasta/nucl defined for sample", sample)
            elif self.params["scope"]=="project":
                if "fasta.nucl" not in self.sample_data:
                    raise AssertionExcept("No fasta/nucl defined for project")
            else:
                raise AssertionExcept("Please supply a scope parameter: either 'sample' or 'project'!")
        else:
            if "scope" not in self.params:
                self.params["scope"]="project"
            elif self.params["scope"]=="sample":
                self.write_warning("It makes no sense to define a sample-scope external reference!")
            elif self.params["scope"]=="project":
                pass
            else:
                self.params["scope"]="project"
        
        
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
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")


        for sample in sample_list:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            use_dir = self.local_start(sample_dir)
            sample_title = sample if sample != "project_data" else self.sample_data["Title"]

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            if "reference" in self.params:
                reference_fasta_file = self.params["reference"]
            else:
                reference_fasta_file = self.sample_data[sample]["fasta.nucl"]

            reference_name = "{dir}{ref_name}_rsem_ref".format(dir=use_dir,ref_name=sample)

            # Get constant part of script:
            self.script += self.get_script_const()
            for other_file in ["gtf","gff3","transcript-to-gene-map","allele-to-gene-map","no-polyA-subset"]:
                if "--"+other_file not in self.params["redir_params"]:  # If passed by user, do not include internal
                    if other_file in self.sample_data[sample]:          # If file exists internally
                        self.script += "--{tag} {value} \\\n\t".format(tag=other_file,
                                                                       value=self.sample_data[sample][other_file])
            self.script += "%s \\\n\t"  % reference_fasta_file
            self.script += "%s \n\n"  % reference_name

            self.sample_data[sample]["RSEM_index"] = "{dir}{ref_name}_rsem_ref".format(dir=sample_dir, ref_name=sample)
            self.sample_data[sample]["RSEM_fasta"] = reference_fasta_file
            if "--star" in self.params["redir_params"]:
                self.sample_data[sample]["STAR.index"] = "{dir}".format(dir=sample_dir)
                self.sample_data[sample]["STAR.fasta"] = reference_fasta_file

            # TODO: Check the exact part of the name to include in the reference!!
            if "--bowtie" in self.params["redir_params"]:
                raise AssertionExcept("RSEM prep with bowtie is not yet defined. Please help with this!")
                # self.sample_data[sample]["bowtie1_index"] = "{dir}".format(dir=sample_dir)
                # self.sample_data[sample]["bowtie1_fasta"] = reference_fasta_file

            if "--bowtie2" in self.params["redir_params"]:
                raise AssertionExcept("RSEM prep with bowtie2 is not yet defined. Please help with this!")
                # self.sample_data[sample]["bowtie2_index"] = "{dir}".format(dir=sample_dir)
                # self.sample_data[sample]["bowtie2_fasta"] = reference_fasta_file

            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
     
        # else:    #if self.params["scope"] == "project":
        #
        #
        #     # Name of specific script:
        #     self.spec_script_name = self.set_spec_script_name()
        #     self.script = ""
        #
        #
        #     # This line should be left before every new script. It sees to local issues.
        #     # Use the dir it returns as the base_dir for this step.
        #     use_dir = self.local_start(self.base_dir)
        #
        #     if "reference" in self.params:
        #         reference_fasta_file = self.params["reference"]
        #     else:
        #         reference_fasta_file = self.sample_data["project_data"]["fasta.nucl"]
        #
        #     reference_name = "{dir}{ref_name}_rsem_ref".format(dir=use_dir,ref_name=self.sample_data["Title"])
        #
        #     # Get constant part of script:
        #     self.script += self.get_script_const()
        #     for other_file in ["gtf","gff3","transcript-to-gene-map","allele-to-gene-map","no-polyA-subset"]:
        #         if "--%s"%other_file not in self.params["redir_params"]:  # If passed by user, do not include internal
        #             if other_file in self.sample_data:          # If file exists internally
        #                 self.script += "--{tag} {value} \\\n\t".format(tag=other_file, value=self.sample_data["project_data"][other_file])
        #     self.script += "%s \\\n\t"  % reference_fasta_file
        #     self.script += "%s \n\n"  % reference_name
        #
        #
        #     self.sample_data["project_data"]["RSEM_index"] = "{dir}{ref_name}_rsem_ref".format(dir=self.base_dir,ref_name=self.sample_data["Title"])
        #     self.sample_data["project_data"]["RSEM_fasta"] = reference_fasta_file
        #
        #     if "--star" in self.params["redir_params"]:
        #         self.sample_data["project_data"]["STAR.index"] = "{dir}".format(dir=self.base_dir)
        #         self.sample_data["project_data"]["STAR.fasta"] = reference_fasta_file
        #
        #     self.local_finish(use_dir,self.base_dir)
        #     self.create_low_level_script()