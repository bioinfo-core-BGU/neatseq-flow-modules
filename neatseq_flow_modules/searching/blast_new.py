# -*- coding: UTF-8 -*-
""" 
``blast_new`` - Alternative implementation
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for executing BLAST of any type on a nucleotide or protein fasta file.
The search can be either on a sample fasta or on a project-wide fasta.
It can use the fasta as a database or as a query.
If used as a database, you must call the makeblastdb module prior to this step.

both ``-query`` and ``-db`` must be passed in the redirected parameters. They should be set to one of the following values:

* ``sample`` - The ``-query`` or ``-db`` should be taken from the sample scope

* ``project`` - The ``-query`` or ``-db`` should be taken from the project scope
  
* **A path** - A path to a ``fasta`` file or ``makeblastdb`` database to use in the parameter as-is.

The type of fasta and database to use are set with the ``querytype`` and ``dbtype`` parameters, respectively. ``dbtype`` must be set if ``-db`` is set to ``sample`` or ``project``. **``querytype`` must be set regardless.** It will determine the type of blast report (*i.e.* whether it will be stored in ``blast.nucl`` or ``blast.prot``)


    
Requires:
~~~~~~~~~~~~~

* fasta files in one of the following slots for sample-wise blast:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

* or fasta files in one of the following slots for project-wise blast:

    * ``sample_data["fasta.nucl"]``
    * ``sample_data["fasta.prot"]``
    
* or a ``makeblastdb`` index in one of the following slots:

    * When ``-db`` is set to 'project'
    
        * ``sample_data["blastdb.nucl"|"blastdb.prot"]``

    * When ``-db`` is set to 'sample'

        * ``sample_data[<sample>]["blastdb.nucl"|"blastdb.prot"]``

    
    
Output:
~~~~~~~~~~~~~

* puts BLAST output files in the following slots for sample-wise blast:

    * ``sample_data[<sample>]["blast.nucl"|"blast.prot"]``
    * ``sample_data[<sample>]["blast"]``

* puts fasta output files in the following slots for project-wise blast:
    
    * ``sample_data["blast.nucl"|"blast.prot"]``
    * ``sample_data["blast"]``


Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10
    
    "dbtype", "nucl|prot", "Helps the module decide which blastdb to use."
    "querytype", "nucl|prot", "Helps the module decide which fasta file to use."
    "redirects: -query | -db", "sample|project|<Path to fasta or BLAST index>", "Redirected compulsory parameters. Set to ``sample`` for sample-scope files, to ``project`` for project-scope files, or to a path for an external file."

.. Note:: You can't set both ``-db`` and ``-query`` to external files. One of them at least has to be ``sample`` or ``project``.


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


External query, project-wise *nucl*-type database (must be proceeded by ``makeblastdb`` module)::


    tbl_blst_int:
        module: blast_new
        base: mkblst1
        script_path: /path/to/blastn
        dbtype: nucl
        redirects:
            -query: /path/to/query.fasta
            -db:  project
            -evalue: 0.0001
            -max_target_seqs: 5
            -num_of_proc: 20
            -num_threads: 20


Sample specific *prot*-type fasta, external database::

    tbl_blst_ext:
        module: blast_new
        base: prokka1
        script_path: /path/to/blastp
        querytype: prot
        redirects:
            -db: /path/to/blasdb.ind
            -query: sample
            -evalue: 0.0001
            

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. **Gapped BLAST and PSI-BLAST: a new generation of protein database search programs**. *Nucleic acids research*, 25(17), pp.3389-3402.

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_blast_new(Step):
    
    auto_redirs = "-db -query".split(" ")

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".blast.out"

        # # Check that either -db or -query (not both) are set in redir_params:
        # if "-db" not in self.params["redir_params"]:
            # raise AssertionExcept("You must supply a '-db' redirects parameter\n")
        # if "-query" not in self.params["redir_params"]:
            # raise AssertionExcept("You must supply a '-query' redirects parameter\n")
        # Check that either db or query (not both) are set:
        if "db" not in self.params:
            raise AssertionExcept("You must supply a 'db' parameter\n")
        if "query" not in self.params:
            raise AssertionExcept("You must supply a 'query' parameter\n")
        
        
        
        if self.params["db"] in ["sample","project"]:
            if "dbtype" not in self.params:
                raise AssertionExcept("No 'dbtype' passed. Please specify the db type to use.")
            elif self.params["dbtype"] not in ["prot","nucl"]:
                raise AssertionExcept("'dbtype' must be either 'prot' or 'nucl'.")
            else:
                pass

        if "querytype" not in self.params:
            raise AssertionExcept("No 'querytype' passed. Please specify the query type to use. (If external query, define it's type with 'querytype')")
        elif self.params["querytype"] not in ["prot","nucl"]:
            raise AssertionExcept("'querytype' must be either 'prot' or 'nucl'.")
        else:
            pass

        if self.params["query"] not in ["sample","project"] and self.params["db"] not in ["sample","project"]:
            raise AssertionExcept("At least one of 'query' and 'db' must be 'sample' or 'project'. You can't pass two paths.")

            

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        
        # Store output format
        outfmt_txt = ["pairwise", "query-anchored showing identities", "query-anchored no identities", "flat query-anchored, show identities", "flat query-anchored, no identities", "XML Blast output", "tabular", "tabular with comment lines", "Text ASN.1", "Binary ASN.1", "Comma-separated values", "BLAST archive format (ASN.1)"]
        if "-outfmt" in self.params["redir_params"]:
            if isinstance(self.params["redir_params"]["-outfmt"], str):
                num_type = re.search("[\'\"]*(\d+)",self.params["redir_params"]["-outfmt"])
                try:
                    
                    blast_type_num = int(num_type.groups()[0])
                    if not 0 <= blast_type_num <= 18:
                        raise AssertionExcept("-outfmt must be between 0 and 18\n")
                    self.outfmt = blast_type_num
                except:
                    raise
            elif isinstance(self.params["redir_params"]["-outfmt"], int):
                if not 0 < self.params["redir_params"]["-outfmt"] < 18:
                    raise AssertionExcept("-outfmt must be between 0 and 18\n")
                self.outfmt = self.params["redir_params"]["-outfmt"]
            else:
                raise AssertionExcept("Unrecognized -outfmt format")
        else:
            self.outfmt = 0
        
        self.outfmt_txt = outfmt_txt[self.outfmt]

        
        if self.params["query"] == "sample" or self.params["db"] == "sample":
            self.params["scope"] = "sample"
            self.step_sample_initiation_bysample()
        else:
            self.params["scope"] = "project"
            self.step_sample_initiation_byproject()

        # Moving the redirected values of -query and -db to self.
        self.db = self.params["db"]
        self.query = self.params["query"]
        # # Removing from redirects. Get treated individually.
        # del self.params["redir_params"]["-db"]
        # del self.params["redir_params"]["-query"]
        
        
        
    def step_sample_initiation_bysample(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for sample-level BLAST
        """
        
            
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            if self.params["query"] == "sample":
                if "fasta." + self.params["querytype"] not in self.sample_data[sample]:
                    raise AssertionExcept("No fasta of type %s in sample" % self.params["querytype"], sample)
            if self.params["query"] == "project":
                if "fasta." + self.params["querytype"] not in self.sample_data:
                    raise AssertionExcept("No fasta of type %s in project" % self.params["querytype"])
            if self.params["db"] == "sample":
                if "blastdb." + self.params["dbtype"] not in self.sample_data[sample]:
                    raise AssertionExcept("No blastdb of type %s in sample. Did you run makeblastdb module with sample scope?" % self.params["dbtype"], sample)            
            if self.params["db"] == "project":
                if "blastdb." + self.params["dbtype"] not in self.sample_data:
                    raise AssertionExcept("No blastdb of type %s in project. Did you run makeblastdb module with project scope?" % self.params["dbtype"])
            
       
            self.sample_data[sample]["blast.outfmt"]     = self.outfmt
            self.sample_data[sample]["blast.outfmt.txt"] = self.outfmt_txt

        
        
        

    def step_sample_initiation_byproject(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for project-level BLAST
        """
        

        if self.params["query"] == "project":
            if "fasta." + self.params["querytype"] not in self.sample_data:
                raise AssertionExcept("No fasta of type %s in project" % self.params["querytype"])
        if self.params["db"] == "project":
            if "blastdb." + self.params["dbtype"] not in self.sample_data:
                raise AssertionExcept("No blastdb of type %s in project. Did you run makeblastdb module with project scope?" % self.params["dbtype"])
                   
        self.sample_data["blast.outfmt"]     = self.outfmt
        self.sample_data["blast.outfmt.txt"] = self.outfmt_txt

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """

        

          
        if self.params["scope"]=="project":
            pass
        elif self.params["scope"]=="sample":
            self.make_sample_file_index()   # see definition below
        
        

    
    def build_scripts(self):
        """ This is the actual script building function
            
        """
        

          
        if self.params["scope"]=="project":
            self.build_scripts_byproject()
        elif self.params["scope"]=="sample":
            self.build_scripts_bysample()
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")
                
        

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
        """
   
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
                
            # Define output filename 
            output_filename = "".join([use_dir , sample , self.file_tag])

            self.script += self.get_script_const()

            # Adding -db :
            ## If -db scope is 'project':
            if self.db == "project":
                ## If dbtype is specified by user, look for the equivalent db in blastdb.nucl or blastdb.prot (both set by makeblastdb module)
                if "dbtype" in self.params:
                    try:
                        self.script += "-db %s \\\n\t" % self.sample_data["blastdb." + self.params["dbtype"]]
                    except KeyError:
                        raise AssertionExcept("No blastdb of type %s exists" % self.params["dbtype"])
                ## If dbtype is NOT specified by user, use default blastdb set by makeblastdb. Let the user beware...
                else:
                    self.script += "-db %s \\\n\t" % self.sample_data["blastdb"]
            ## Same as above but for -db in sample scope.
            elif self.db == "sample":
                if "dbtype" in self.params:
                    try:
                        self.script += "-db %s \\\n\t" % self.sample_data[sample]["blastdb." + self.params["dbtype"]]
                    except KeyError:
                        raise AssertionExcept("No blastdb of type %s exists" % self.params["dbtype"], sample)
                else:
                    self.script += "-db %s \\\n\t" % self.sample_data[sample]["blastdb"]
            else: # -db is a user-defined path:
                self.script += "-db %s \\\n\t" % self.db
                    
            # Adding -query :
            ## Same as for -db. See documentation above
            if self.query == "project":
                if "querytype" in self.params:
                    try:
                        self.script += "-query %s \\\n\t" % self.sample_data["fasta." + self.params["querytype"]]
                    except KeyError:
                        raise AssertionExcept("No fasta of type %s exists" % self.params["querytype"])
                else:
                    raise AssertionExcept("You must specify querytype")
            elif self.query == "sample":
                if "querytype" in self.params:
                    try:
                        self.script += "-query %s \\\n\t" % self.sample_data[sample]["fasta." + self.params["querytype"]]
                    except KeyError:
                        raise AssertionExcept("No blastdb of type %s exists" % self.params["querytype"], sample)
                else:
                    raise AssertionExcept("You must specify querytype", sample)
            else: # Path
                self.script += "-query %s \\\n\t" % self.query
                
            self.script += "-out %s\n\n" % output_filename
            
            # Store BLAST result file:
            self.sample_data[sample]["blast"] = (sample_dir + os.path.basename(output_filename))
            self.sample_data[sample]["blast." + self.params["querytype"]] = self.sample_data[sample]["blast"]
            self.stamp_file(self.sample_data[sample]["blast"])
            
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            
            self.create_low_level_script()
                    
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """


        
        
        # Each iteration must define the following class variables:
        # spec_script_name
        # script
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
                
                
        # Define output filename 
        output_filename = "".join([use_dir , self.sample_data["Title"] , self.file_tag])

        self.script += self.get_script_const()
        # Adding -db :
        if self.db == "project":
            if "dbtype" in self.params:
                try:
                    self.script += "-db %s \\\n\t" % self.sample_data["blastdb." + self.params["dbtype"]]
                except KeyError:
                    raise AssertionExcept("No blastdb of type %s exists" % self.params["dbtype"])
            else:
                self.script += "-db %s \\\n\t" % self.sample_data["blastdb"]
        else: # Path
            self.script += "-db %s \\\n\t" % self.db
                
        # Adding -query :
        if self.query == "project":
            if "querytype" in self.params:
                try:
                    self.script += "-query %s \\\n\t" % self.sample_data["fasta." + self.params["querytype"]]
                except KeyError:
                    raise AssertionExcept("No fasta of type %s exists" % self.params["querytype"])
            else:
                raise AssertionExcept("You must specify querytype")
        else: # Path
            self.script += "-query %s \\\n\t" % self.query
            
        self.script += "-out %s\n\n" % output_filename
            
        # Store BLAST result file:
        self.sample_data["blast"] = (self.base_dir + os.path.basename(output_filename))
        self.sample_data["blast." + self.params["querytype"]] = self.sample_data["blast"]
        self.stamp_file(self.sample_data["blast"])



        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
                

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names.
            This can be used by scripts called by create_spec_wrapping_up_script() to summarize the BLAST outputs.
        """
        
        with open(self.base_dir + "BLAST_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tBLAST_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["blast." + self.fasta2use]))
                
        self.sample_data["BLAST_files_index"] = self.base_dir + "BLAST_files_index.txt"
        
  
        
