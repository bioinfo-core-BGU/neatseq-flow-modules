# -*- coding: UTF-8 -*-
""" 
``RnammerTranscriptome``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for executing ``RnammerTranscriptome.pl`` script from Trinotate to use RNAMMER for searching for rRNA
sequernces.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* If ``scope = sample``, ``fasta`` files in the following slot:
    
    * ``sample_data[<sample>]["fasta.nucl"]``

* If ``scope = project``, ``fasta`` files in the following slot:
    
    * ``sample_data["fasta.nucl"]``

    
Output:
~~~~~~~~~~~~~

* puts the resulting gff file in the following slots::
        
    * for ``scope = sample``:
    
        * ``sample_data[<sample>]["gff"]``
        * ``sample_data[<sample>]["rnammer"]``
    
    * for ``scope = project``:
    
        * ``sample_data["gff"]``
        * ``sample_data["rnammer"]``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Create one assembly for all samples or one assembly per sample."

.. Attention:: Please set ``--path_to_rnammer`` in redirects to the path to RNAMMER executable.
    See the Trinotate documentation for installation instructions (it's a bit tricky

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_rnammer_highExpr:
        module:             RnammerTranscriptome
        base:               Split_Fasta
        script_path:        {Vars.paths.RnammerTranscriptome}
        scope:              sample
        redirects:
            --path_to_rnammer:  {Vars.paths.rnammer}

            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RNammer: consistent annotation of rRNA genes in genomic sequences Lagesen K, Hallin PF, Rodland E, Staerfeldt HH, Rognes T Ussery DW Nucleic Acids Res. 2007 Apr 22.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_RnammerTranscriptome(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        if "path_to_rnammer" in self.params:
            self.params["redir_params"]["--path_to_rnammer"] = self.params["path_to_rnammer"]
        if "--path_to_rnammer" not in self.params["redir_params"]:
            raise AssertionExcept("Please supply '--path_to_rnammer' via redirects")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                if "fasta.nucl" not in self.sample_data:
                    raise AssertionExcept("Project does not have a 'fasta.nucl' file")

            elif self.params["scope"] == "sample":
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    if "fasta.nucl" not in self.sample_data[sample]:
                        raise AssertionExcept("Sample does not have a 'fasta.nucl' file.", sample)
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")

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
        self.spec_script_name = self.set_spec_script_name()

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        output_basename = "{fasta}.rnammer.gff".format(fasta=os.path.basename(self.sample_data["fasta.nucl"]))

        self.script = """\
cd {dir}

{const}
    --transcriptome {fasta_nucl}
    
cd -
""".format(dir=use_dir,
           const=self.get_script_const().rstrip("\n\t"),
           fasta_nucl=self.sample_data["fasta.nucl"])

        # Store results to fasta and assembly slots:
        self.sample_data["rnammer"] = "{dir}{name}".format(dir=self.base_dir, name=output_basename)
        self.sample_data["gff"] = self.sample_data["rnammer"]

        self.stamp_file(self.sample_data["rnammer"])
        
        # Move all files from temporary local dir to permanent base_dir
        # Sees to copying local files to final destination (and other stuff)
        self.local_finish(use_dir,self.base_dir)

        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            output_basename = "{fasta}.rnammer.gff".\
                format(fasta=os.path.basename(self.sample_data[sample]["fasta.nucl"]))

            self.script = """\
cd {dir}

{const}
    --transcriptome {fasta_nucl}

cd -
""".format(dir=use_dir,
           const=self.get_script_const().rstrip("\n\t"),
           fasta_nucl=self.sample_data[sample]["fasta.nucl"])

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["rnammer"] = "{dir}{name}".format(dir=sample_dir,name=output_basename)
            self.sample_data[sample]["gff"] = self.sample_data[sample]["rnammer"]

            self.stamp_file(self.sample_data[sample]["rnammer"])

            # Wrapping up function. Leave these lines at the end of every iteration:
            # Sees to copying local files to final destination (and other stuff)
            self.local_finish(use_dir,sample_dir)

            self.create_low_level_script()
