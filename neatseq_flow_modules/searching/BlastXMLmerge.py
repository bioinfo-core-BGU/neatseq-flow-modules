# -*- coding: UTF-8 -*-
""" 
``BlastXMLmerge``
-----------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for running ``BlastXMLmerge.py``, available from github at https://github.com/jhbadger/scripts/blob/master/BlastXMLmerge.py

The program merges sample-wise XML blast reports into a single project-wide XML BLAST report.

Can be used together with the ``split_fasta`` module to parallelize BLAST searches: If you have a project-wide fasta file which you want to BLAST, split it into samples with ``split_fasta``, run BLAST on sample scope, and then merge the individual reports with this module.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* XML BLAST result files in:

    * ``sample_data[<sample>]["blast.nucl|blast.prot"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the merged report in:

    * ``self.sample_data["project_data"]["blast.nucl|blast.prot"]``            if ``scope = project``

    

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "blast2use", "``nucl|prot``", "If both nucl and prot BLAST reports exist, you have to specify which one to use with this parameter. If unspecified, will merge both."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    merge_xml_blast:
        module:         BlastXMLmerge
        base:           blst_xml
        script_path:    {Vars.paths.BlastXMLmerge.py}
        blast2use:      nucl

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_BlastXMLmerge(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".blast.parsed"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            if not "blast" in self.sample_data[sample]:
                raise AssertionExcept("There are no BLAST results.\n" , sample)
            if not "blast.outfmt" in self.sample_data[sample]:
                self.write_warning("No 'blast.outfmt.txt' exists for sample. Assuming is XML!")
            if self.sample_data[sample]["blast.outfmt"] not in [5,16]:
                raise AssertionExcept("BLAST report does not seem to be in XML (5) or single-file XML2 (16) format")
            
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    def build_scripts(self):
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
            
            
        # Define output filename 
        output_filename = "".join([self.sample_data["Title"] , self.file_tag])

        # Define query and db files:
        # If db is defined by user, set the query to the correct 'fasta2use'
        # If both nucl and prot appear in blast results
        if "blast.nucl" in self.sample_data and "blast.prot" in self.sample_data:
            if "blast2use" in self.params.keys() and self.params["blast2use"] in ("nucl","prot"):
                blast2use = "blast.%s" % self.params["blast2use"]
                # self.script += "--blast %s \\\n\t" % self.sample_data[sample]["blast"][fasta2use]
            else:
                self.write_warning("Project has both 'nucl' and 'prot' blast results. Using the last one produced.")
                blast2use = "blast"
        elif "blast.nucl" in self.sample_data:
            blast2use = "blast.nucl"
        elif "blast.prot" in self.sample_data:
            blast2use = "blast.prot"
        else:
            raise AssertionExcept("No BLAST Results defined\n")

        
        self.script += self.get_script_const()
        self.script += "{dir}{file}.blast.xml \\\n\t".format(dir = use_dir, file=self.sample_data["Title"])
        for sample in self.sample_data["samples"]:
            self.script += "%s \\\n\t" % self.sample_data[sample][blast2use]
        
        self.script.rstrip("\\\n\t")
        self.script += "\n\n"
        
        # Store BLAST result file:
        self.sample_data["project_data"]["blast"] = "{dir}{file}.blast.xml".format(dir = self.base_dir, file=self.sample_data["Title"])
        if blast2use != "blast":
            self.sample_data["project_data"][blast2use] = self.sample_data["project_data"]["blast"]
        self.stamp_file(self.sample_data["project_data"]["blast"])

        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
                        
                
                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        pass