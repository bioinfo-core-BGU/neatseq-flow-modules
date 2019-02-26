# -*- coding: UTF-8 -*-
""" 
``hmmscan`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for searching a fasta file with hmmscan.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* If ``scope = sample``, ``fasta`` files in one of the following slots: 
    
    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

* If ``scope = project``, ``fasta`` files in one of the following slots: 
    
    * ``sample_data["fasta.nucl"]``
    * ``sample_data["fasta.prot"]``

    
Output:
~~~~~~~~~~~~~

* puts ``hmmscan`` output files in the following slots:
        
    * for ``scope = sample`` (depending on ``type`` passed):
    
        * ``sample_data[<sample>]["hmmscan.nucl"]``
        * ``sample_data[<sample>]["hmmscan.prot"]``
    
    * for ``scope = project`` (depending on ``type`` passed):
    
        * ``sample_data["hmmscan.nucl"]``
        * ``sample_data["hmmscan.prot"]``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Create one assembly for all samples or one assembly per sample."
    "type", "", "Use a `prot` or `nucl` fasta file for the search."
    "output_type", "tblout|domtblout|pfamtblout", "tblout: parseable table of per-sequence hits to file, domtblout: parseable table of per-domain hits to file, pfamtblout: table of hits and domains in Pfam format"
    "hmmdb","","A path to the hmmdb to search against."
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_hmmscan1_highExpr:
        module:             hmmscan
        base:               trino_Transdecode_highExpr
        script_path:        {Vars.paths.hmmscan}
        scope:              sample
        type:               prot
        output_type:        domtblout 
        hmmdb:              {Vars.databases.trinotate.pfam}
        qsub_params:
            -pe:            shared 10
        redirects:
            --cpu:          1

            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finn, Robert D., Jody Clements, and Sean R. Eddy. "HMMER web server: interactive sequence similarity searching." Nucleic acids research 39.suppl_2 (2011): W29-W37.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_hmmscan(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        
        
        if "type" not in self.params:
            raise AssertionExcept("Please specify the fasta type to use: type = nucl or prot")
        if self.params["type"] not in ["nucl","prot"]:
            raise AssertionExcept("'type' must be 'nucl' or 'prot'")
            
        if "output_type" not in self.params:
            raise AssertionExcept("""
Please specify the output_type to use: 
* tblout     : save parseable table of per-sequence hits 
* domtblout  : save parseable table of per-domain hits 
* pfamtblout : save table of hits and domains in Pfam format 
""")
        if self.params["output_type"] not in ["tblout","domtblout","pfamtblout"]:
            raise AssertionExcept("""
'output_type' must be one of the following: 
* tblout     : save parseable table of per-sequence hits 
* domtblout  : save parseable table of per-domain hits 
* pfamtblout : save table of hits and domains in Pfam format 
""")            

        if "hmmdb" not in self.params:
            raise AssertionExcept("Please specify the hmmdb to use!")

        
        
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
            if "fasta.{type}".format(type=self.params["type"]) not in self.sample_data[sample]:
                raise AssertionExcept("No 'fasta.{type}' defined".format(type=self.params["type"]), sample)

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):



        if self.params["scope"]=="project":
            sample_list = ["project_data"]
        elif self.params["scope"]=="sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            output_basename = "{title}.hmmscan.{type}".format(title=sample,
                                                              type=self.params["output_type"])

            self.script = self.get_script_const()
            self.script += "--{output_type} {outfile} \\\n\t".format(output_type=self.params["output_type"],
                                                                     outfile=use_dir+output_basename)
            self.script += "{hmmdb} \\\n\t".format(hmmdb = self.params["hmmdb"])
            self.script += "{seqfile} \\\n\t".format(seqfile = self.sample_data[sample]["fasta.{type}".format(type=self.params["type"])])
            self.script += "> {outfile}.log \n\n".format(outfile = use_dir+output_basename)

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["hmmscan.%s" % self.params["type"]] = "%s%s" % (sample_dir, output_basename)
            self.sample_data[sample]["hmmscan.%s.log" % self.params["type"]] = "%s%s.log" % (sample_dir, output_basename)

            self.stamp_file(self.sample_data[sample]["hmmscan.%s" % self.params["type"]])
            self.stamp_file(self.sample_data[sample]["hmmscan.%s.log" % self.params["type"]])

            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()


