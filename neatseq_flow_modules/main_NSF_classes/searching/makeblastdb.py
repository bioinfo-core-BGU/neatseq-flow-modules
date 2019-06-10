# -*- coding: UTF-8 -*-
""" 
``makeblastdb`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Create a blastdb from a fasta file

        
Requires
~~~~~~~~
        
* fastq files in the following slots:

    * ``sample_data[<sample>]["fasta.nucl"|"fasta.prot"]``

* Or (if 'projectBLAST' is set)
    * ``sample_data["fasta.nucl"|"fasta.prot"]``

Output
~~~~~~
    
* A BLAST database in the following slots:

    * ``sample_data[<sample>]["blastdb"]``
    * ``sample_data[<sample>]["blastdb.nucl"|"blastdb.prot"]``
    * ``sample_data[<sample>]["blastdb.nucl.log"|"blastdb.prot.log"]``

* Or (if 'projectBLAST' is set):

    * ``sample_data["blastdb"]``
    * ``sample_data["blastdb.nucl"|"blastdb.prot"]``
    * ``sample_data["blastdb.nucl.log"|"blastdb.prot.log"]``
    
Parameters that can be set:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide or sample fasta slot should be used"
    "-dbtype", "nucl/prot", "This is a compulsory redirected parameter.Helps the module decide which fasta file to use."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    mkblst1:
        module: makeblastdb
        base: trinity1
        script_path: /path/to/bin/makeblastdb
        redirects:
            -dbtype: nucl
        scope: project

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. **Gapped BLAST and PSI-BLAST: a new generation of protein database search programs**. *Nucleic acids research*, 25(17), pp.3389-3402.
    
"""
    

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_makeblastdb(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "makeblastdb.out"

        # Checking this once and then applying to each sample:
        if not "-dbtype" in self.params["redir_params"]:
            raise AssertionExcept("You must define a -dbtype parameter\n")
        self.dbtype = self.params["redir_params"]["-dbtype"]
        if not self.dbtype in ["nucl","prot"]:
            raise AssertionExcept("-dbtype must be either nucl or prot\n")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if not "scope" in self.params:
            raise AssertionExcept("No 'scope' specified.")

        if self.params["scope"]=="project":
            sample_list = ["project_data"]
        elif self.params["scope"]=="sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        # Creating holder for output:
        for sample in sample_list:      # Getting list of samples out of samples_hash
            # Make sure a file exists in the sample equivalent to dbtype:
            try:
                # In version 1.0.2, nucl and prot slots have been renamed to fasta.nucl and fasta.prot
                self.sample_data[sample]["fasta." + self.dbtype]
            except KeyError:
                raise AssertionExcept("No file exists in sample for specified -dbtype (%s)\n" % self.dbtype, sample)

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        if self.params["scope"]=="project":
            sample_list = ["project_data"]
        elif self.params["scope"]=="sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:
        
            sample_title = sample if sample != "project_data" else self.sample_data["Title"]

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

         
            # Define output filename 
            output_filename = ".".join([sample_title, self.name, self.file_tag])
            blastdb_title = os.path.basename(output_filename)

            self.script += self.get_script_const()
            self.script += "-out {dir}{fn} \\\n\t".format(dir= use_dir, fn=output_filename)
            self.script += "-in %s \\\n\t" % self.sample_data[sample]["fasta." + self.dbtype]
            self.script += "-title %s \\\n\t" % blastdb_title
            self.script += "-logfile {dir}{outname}.log \n\n".format(dir=use_dir,outname=output_filename)

            self.sample_data[sample]["blastdb." + self.dbtype] = (sample_dir + output_filename)
            self.sample_data[sample]["blastdb." + self.dbtype + ".log"] = "{dir}{fn}.log".format(dir= sample_dir,
                                                                                                 fn=output_filename)

            self.sample_data[sample]["blastdb"] = self.sample_data[sample]["blastdb." + self.dbtype]
            
            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()

