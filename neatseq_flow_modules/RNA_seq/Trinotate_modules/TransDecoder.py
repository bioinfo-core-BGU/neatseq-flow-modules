# -*- coding: UTF-8 -*-
""" 
``TransDecoder`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``TransDecoder`` on a transcripts file.

.. Note:: Tested on TransDecoder version 5.5.0.. The main difference being that in this version an output directory can be specified in the command line.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
``fasta`` files in at least one of the following slots:
    
    * ``sample_data[<sample>]["fasta.nucl"]``  (if ``scope = sample``)
    * ``sample_data["fasta.nucl"]``  (if ``scope = project``)

    
Output:
~~~~~~~~~~~~~

* If ``scope = project``:

    * Protein fasta in ``self.sample_data["project_data"]["fasta.prot"]``
    * Gene fasta in ``self.sample_data["project_data"]["fasta.nucl"]``
    * Original transcripts in ``self.sample_data["project_data"]["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data["project_data"]["gff3"]``

* If ``scope = sample``:

    * Protein fasta in ``self.sample_data[<sample>]["fasta.prot"]``
    * Gene fasta in ``self.sample_data[<sample>]["fasta.nucl"]``
    * Original transcripts in ``self.sample_data[<sample>]["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data[<sample>]["gff3"]``


                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Determine weather to use sample or project transcripts file."
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_Transdecode_highExpr:
        module:             TransDecoder
        base:               Split_Fasta
        script_path:        {Vars.paths.TransDecoder}
        scope:              sample
        
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_TransDecoder(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".Trinity.fasta"
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        if self.params["scope"]=="project":
            sample_list = ["project_data"]
        elif self.params["scope"]=="sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:
            if "fasta.nucl" not in self.sample_data[sample]:
                raise AssertionExcept("No 'fasta.nucl' defined!", sample)

            if "Predict" in self.params or re.search("Predict", self.params["script_path"]):
                # Adding directory from previous run:
                if "transdecoder.dir" not in self.sample_data[sample]:
                    raise AssertionExcept("Please include a 'LongOrf' TransDecoder step before the 'Predict' step.\n"
                                          "You can include 'blatsp' and 'hmmscan' steps in between to make it effective")

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):
    
        # if self.params["scope"] == "project":
        #     self.build_scripts_project()
        # else:
        #     self.build_scripts_sample()

        if self.params["scope"]=="project":
            sample_list = ["project_data"]
        elif self.params["scope"]=="sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:
            sample_title = self.set_sample_name(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)


            # self.script = "cd {usedir}\n\n".format(usedir=use_dir)
            self.script += self.get_script_const()

            if "Predict" in self.params or re.search("Predict",self.params["script_path"]):
                # Adding directory from previous run:
                self.script += "-O {outdir} \\\n\t".format(outdir=self.sample_data[sample]["transdecoder.dir"])
                # Adding hmm and blastp searches if exist:
                if not "hmmscan.prot" in self.sample_data[sample] and not "blast.prot" in self.sample_data[sample]:
                    raise AssertionExcept("Cannot predict without either hmmscan of pfam or blastp.")
                else:
                    if "hmmscan.prot" in self.sample_data[sample]:  # and not "blast.prot" in self.sample_data[sample]:
                        self.script += "--retain_pfam_hits {hmmpfam} \\\n\t".format(
                            hmmpfam=self.sample_data[sample]["hmmscan.prot"])
                    if "blast.prot" in self.sample_data[sample]:
                        self.script += "--retain_blastp_hits {blastp} \\\n\t".format(
                            blastp=self.sample_data[sample]["blast.prot"])

            else:       # LongOrfs
                # If LongOrf, set the TransDecoder directory. Will be stored and used by Predict
                # output_dir_base = "%s.transdecoder_dir" % os.path.basename(self.sample_data[sample]["fasta.nucl"])
                self.script += "-O {dir} \\\n\t".format(dir=use_dir)
                self.sample_data[sample]["transdecoder.dir"] = sample_dir
                # Seeing to gene_trans_map:
                if "--gene_trans_map" in self.params["redir_params"]:
                    if self.params["redir_params"]["--gene_trans_map"]:  # If value was passed
                        self.sample_data[sample]["gene_trans_map"] = self.params["redir_params"]["--gene_trans_map"]
                    elif "gene_trans_map" in self.sample_data[sample]:  # If passed empty, use internal:
                        self.params["redir_params"]["--gene_trans_map"] = self.sample_data[sample]["gene_trans_map"]
                    else:  # If passed empty, but no internal gene_trans_map exists:
                        raise AssertionExcept("Expecting 'gene_trans_map' but none found.\n", sample)

            self.script += "-t {fasta_nucl}\n\n".format(fasta_nucl=self.sample_data[sample]["fasta.nucl"])

            # self.script += "cd {home_dir}\n\n".format(home_dir=self.pipe_data["home_dir"])

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["fasta.prot"] = os.path.join(sample_dir, "longest_orfs.pep")
            self.sample_data[sample]["transcripts.fasta.nucl"] = self.sample_data[sample]["fasta.nucl"]
            self.sample_data[sample]["fasta.nucl"] = os.path.join(sample_dir, "longest_orfs.cds")
            self.sample_data[sample]["gff3"] = os.path.join(sample_dir, "longest_orfs.gff3")

            self.stamp_file(self.sample_data[sample]["fasta.prot"])
            self.stamp_file(self.sample_data[sample]["fasta.nucl"])
            self.stamp_file(self.sample_data[sample]["gff3"])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir, sample_dir)  # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
