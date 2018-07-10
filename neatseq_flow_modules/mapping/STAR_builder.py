# -*- coding: UTF-8 -*-
""" 
``STAR_builder``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running STAR genome index construction:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.nucl"]``

* If ``--sjdbGTFfile`` is set in redirects, but left empty, will expect to find a ``GTF`` file here:

    * ``sample_data["gtf"]``  if ``scope`` = "project"
    * ``sample_data[<sample>]["gtf"]``   if ``scope`` = `"sample"

* If ``--sjdbFileChrStartEnd`` is set in redirects, but left empty, will expect to find an SJ file here:

    * ``sample_data["SJ.out.tab"]``  if ``scope`` = "project"
    * ``sample_data[<sample>]["SJ.out.tab"]``   if ``scope`` = "sample"


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Puts output index files in one of the following slot:

    * ``self.sample_data[<sample>]["STAR.index"]``
    * ``self.sample_data["STAR.index"]``


Puts the reference fasta file in one of the following slot:

        * ``self.sample_data[<sample>]["STAR.fasta"]``
    * ``self.sample_data["STAR.fasta"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "Not used"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::


    STAR_bld_ind:
        module:             STAR_builder
        base:               trinity1
        script_path:        /path/to/STAR
        scope:              project
        qsub_params:
            queue:          star.q
        redirects:
            --genomeSAindexNbases:  12
            --genomeChrBinNbits:    10

    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*, 29(1), pp.15-21.

"""


import os, re
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_STAR_builder(Step):

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        self.auto_redirs = ["--runMode", "--genomeDir", "--genomeFastaFiles"]

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # By default, do not use internal GTF and SJ files, except when requested in param file
        # See below
        self.use_internal_GTF_file = False
        self.use_internal_SJ_file = False
        if "scope" not in self.params.keys():
            # Try guessing scope:
            try:  # Does a nucl fasta exist for project?
                self.sample_data["fasta.nucl"]
            except KeyError:
                self.params["scope"] = "sample"
            else:
                self.params["scope"] = "project"
        else:
            # Check scope is legitimate
            if not self.params["scope"] in ["project","sample"]:
                raise AssertionExcept("Scope must be either 'sample' or 'project'\n")

        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a nucl fasta defined. Can't build index\n", sample)
                # Checking the the splice junctions file exists, if requested:
                # 1. --sjdbFileChrStartEnd exists in params
                # 2. --sjdbFileChrStartEnd is empty
                if "--sjdbFileChrStartEnd" in self.params["redir_params"] and \
                    not self.params["redir_params"]["--sjdbFileChrStartEnd"]:
                    if "SJ.out.tab" in self.sample_data[sample]:
                        self.use_internal_SJ_file = True
                    else:
                        raise AssertionExcept("You passed an empty '--sjdbFileChrStartEnd' but sample does not have "
                                              "a 'SJ.out.tab' file defined. "
                                              "Can't build index\n", sample)
                if "--sjdbGTFfile" in self.params["redir_params"] and \
                    not self.params["redir_params"]["--sjdbGTFfile"]:
                    if "gtf" in self.sample_data[sample]:
                        self.use_internal_GTF_file = True
                    else:
                        raise AssertionExcept("You passed an empty '--sjdbGTFfile' but sample does not have "
                                              "a 'gtf' file defined. "
                                              "Can't build index\n", sample)
        else:
            try:
                self.sample_data["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("Project does not have a nucl fasta defined. Can't build index\n")
            # Checking the the splice junctions file exists, if requested:
            # 1. --sjdbFileChrStartEnd exists in params
            # 2. --sjdbFileChrStartEnd is empty
            if "--sjdbFileChrStartEnd" in self.params["redir_params"] and \
                not self.params["redir_params"]["--sjdbFileChrStartEnd"]:
                if "SJ.out.tab" in self.sample_data:
                    self.use_internal_SJ_file = True
                else:
                    raise AssertionExcept("You passed an empty '--sjdbFileChrStartEnd' but project does not have "
                                          "a 'SJ.out.tab' file defined. "
                                          "Can't build index\n")
            if "--sjdbGTFfile" in self.params["redir_params"] and \
                    not self.params["redir_params"]["--sjdbGTFfile"]:
                if "gtf" in self.sample_data:
                    self.use_internal_GTF_file = True
                else:
                    raise AssertionExcept("You passed an empty '--sjdbGTFfile' but project does not have "
                                          "a 'gtf' file defined. "
                                          "Can't build index\n")

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
        if self.params["scope"] == "sample":
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                # Requested internal SJ.out.tab file.
                # Setting in redir_params in each instance to sample SJ.out.tab file
                if self.use_internal_SJ_file:
                    self.params["redir_params"]["--sjdbFileChrStartEnd"] = self.sample_data[sample]["SJ.out.tab"]
                # Requested internal gtf file.
                # Setting in redir_params in each instance to sample gtf
                if self.use_internal_GTF_file:
                    self.params["redir_params"]["--sjdbGTFfile"] = self.sample_data[sample]["gtf"]

                # Get constant part of script:
                self.script += self.get_script_const()
                self.script += "--runMode genomeGenerate \\\n\t"
                self.script += "--genomeDir %s \\\n\t"  % use_dir
                self.script += "--genomeFastaFiles %s \n\n"  % self.sample_data[sample]["fasta.nucl"]

                self.sample_data[sample]["STAR.index"] = sample_dir
                self.sample_data[sample]["STAR.fasta"] = self.sample_data[sample]["fasta.nucl"]

                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)
                self.create_low_level_script()
        else:
            self.spec_script_name = self.set_spec_script_name()

            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Requested internal SJ.out.tab file.
            # Setting in redir_params in each instance to sample SJ.out.tab file
            if self.use_internal_SJ_file:
                self.params["redir_params"]["--sjdbFileChrStartEnd"] = self.sample_data["SJ.out.tab"]
            # Requested internal GTF file.
            # Setting in redir_params in each instance to project gtf file
            if self.use_internal_GTF_file:
                self.params["redir_params"]["--sjdbGTFfile"] = self.sample_data["gtf"]

            # Get constant part of script:
            self.script += self.get_script_const()
            self.script += "--runMode genomeGenerate \\\n\t"
            self.script += "--genomeDir %s \\\n\t"  % use_dir
            self.script += "--genomeFastaFiles %s \n\n"  % self.sample_data["fasta.nucl"]

            self.sample_data["STAR.index"] = self.base_dir
            self.sample_data["STAR.fasta"] = self.sample_data["fasta.nucl"]
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
