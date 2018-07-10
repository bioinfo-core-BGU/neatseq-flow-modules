# -*- coding: UTF-8 -*-
""" 
``Trinotate`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for RNA_seq assembly annotation using `Trinotate`_.

.. _Trinotate: https://trinotate.github.io/

.. Note:: This module will be updated in the future to support uploading of other sources of information such as RNAMMER output. See `Trinotate`_ documentation.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* A transcripts file in 
    * self.sample_data["transcripts.fasta.nucl"],
* A gene to transcript mapping file in: (produced by ``Trinity_gene_to_trans_map`` module)
    * self.sample_data["gene_trans_map"], 
* A protein fasta file (produced by ``TransDecoder``)
    * self.sample_data["fasta.prot"])
* Results of ``blastx`` of protein file against swissprot database:
    * self.sample_data["blast.prot"], 
* Results of ``blastn`` of transcripts file against swissprot database:
    * self.sample_data["blast.nucl"], 
* Results of ``hmmscan`` of protein file against pfam database:
    * self.sample_data["hmmscan.prot"])
    
.. Attention:: If ``scope`` is set to ``sample``, all of the above files should be in the sample scope!

    
Output:
~~~~~~~~~~~~~

* puts Trinotate report file in:

    * ``sample_data[<sample>]["trino.rep"]`` (``scope = sample``)
    * ``sample_data["trino.rep"]`` (``scope = project``)
        

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", ""
    "sqlitedb","","Path to Trinotate sqlitedb"
    "cp_sqlitedb","","Create local copy of the sqlitedb, before loading teh data (recommended)"
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_Trinotate:
        module:             Trinotate
        base:               
                            - trino_blastp_sprot
                            - trino_blastx_sprot
                            - trino_hmmscan1
        script_path:        {Vars.paths.Trinotate}
        scope:              project
        sqlitedb:           {Vars.databases.trinotate.sqlitedb}
        cp_sqlitedb:    
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011. **Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data**. *Nature biotechnology*, 29(7), p.644.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_Trinotate(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        
        
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                if "fasta.nucl" not in self.sample_data:
                    raise AssertionExcept("Project does not have a nucl fasta.")
                if "fasta.prot" not in self.sample_data:
                    raise AssertionExcept("Project does not have a prot fasta.")
                if "hmmscan.prot" not in self.sample_data:
                    raise AssertionExcept("Project does not have a prot hmmscan output file.")

            elif self.params["scope"]=="sample":
                
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    if "fasta.nucl" not in self.sample_data[sample]:
                        raise AssertionExcept("Project does not have a nucl fasta.")
                    if "fasta.prot" not in self.sample_data[sample]:
                        raise AssertionExcept("Project does not have a prot fasta.")
                    if "hmmscan.prot" not in self.sample_data[sample]:
                        raise AssertionExcept("Project does not have a prot hmmscan output file.")
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

        output_basename = "{title}.trino_anno_rep.xls".format(title = self.sample_data["Title"])

        # Setting the sqlitedb to use:
        self.sqlitedb = self.params["sqlitedb"]
         
        # If requested copy, create copy and change active sqlitedb
        if "cp_sqlitedb" in self.params:
            self.script += "cp {source} {dest}\n\n".format(source=self.params["sqlitedb"],
                                                           dest=os.path.join(use_dir,os.path.basename(self.params["sqlitedb"])))
            # Setting sqlitedb to new copy of db:
            self.sqlitedb = os.path.join(use_dir, os.path.basename(self.params["sqlitedb"]))
            
        trino_cmd_sqlite = "{script} \\\n\t{sqlitedb} \\\n\t".format(script=self.params["script_path"],
                                                                     sqlitedb=self.sqlitedb)

        ################################ Step 1. init
        self.script += "### Step 1: init db\n\n"
        self.script += """
{trino_cmd}init \\
\t--gene_trans_map {trans_map} \\
\t--transcript_fasta {trans_fa} \\
\t--transdecoder_pep  {pep_fa}\n\n""".format(trino_cmd = trino_cmd_sqlite,
                                             trans_map = self.sample_data["gene_trans_map"],
                                             trans_fa  = self.sample_data["transcripts.fasta.nucl"],
                                             pep_fa    = self.sample_data["fasta.prot"])
        
        
        ################################ Step 2. Load
        self.script += "### Step 2: Load reports\n\n"
        self.script += """
{trino_cmd}LOAD_swissprot_blastp {blastp} \n
{trino_cmd}LOAD_swissprot_blastx {blastx} \n

""".format(trino_cmd = trino_cmd_sqlite,
           blastp  = self.sample_data["blast.prot"],
           blastx = self.sample_data["blast.nucl"])

        if "hmmscan.prot" in self.sample_data:
            self.script += "{trino_cmd}LOAD_pfam {pfam}\n\n".format(trino_cmd = trino_cmd_sqlite,
                                                                    pfam   = self.sample_data["hmmscan.prot"])

        if "rnammer" in self.sample_data:
            self.script += "{trino_cmd}LOAD_rnammer {rnammer}\n\n".format(trino_cmd=trino_cmd_sqlite,
                                                                    rnammer=self.sample_data["rnammer"])

        ################################ Step 4. Report
        self.script += "### Step 4: Create report\n\n"
        self.script += """
{trino_cmd}report \\
\t{redirects} > {dir}{file} \n\n""".format(trino_cmd = trino_cmd_sqlite,
                                         redirects   = self.get_redir_parameters_script(),
                                         dir         = use_dir, 
                                         file        = output_basename)
        

        # Store results to fasta and assembly slots:
        self.sample_data["trino.rep"] = "%s%s" % (self.base_dir, output_basename)
        
        self.stamp_file(self.sample_data["trino.rep"])

        ################################ Step 5. delete sqlitedb
        # If requested copy, create copy and change active sqlitedb
        if "cp_sqlitedb" in self.params and "rm_sqlitedb" in self.params:
            self.script += "rm -rf {dest}\n\n".format(dest=os.path.join(use_dir,os.path.basename(self.params["sqlitedb"])))

            
        # Move all files from temporary local dir to permanent base_dir
        # Sees to copying local files to final destination (and other stuff)
        self.local_finish(use_dir,self.base_dir)

        self.create_low_level_script()

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
            
            output_basename = "{title}.trino_anno_rep.xls".format(title = sample)

            # Setting the sqlitedb to use:
            self.sqlitedb = self.params["sqlitedb"]
             
            # If requested copy, create copy and change active sqlitedb
            if "cp_sqlitedb" in self.params:
                self.script += "cp {source} {dest}\n\n".format(source=self.sqlitedb,
                                                               dest=os.path.join(use_dir,os.path.basename(self.params["sqlitedb"])))
                # Setting sqlitedb to new copy of db:
                
                self.sqlitedb = os.path.join(use_dir,os.path.basename(self.params["sqlitedb"]))
                
            trino_cmd_sqlite = "{script} \\\n\t{sqlitedb} \\\n\t".format(script = self.params["script_path"],
                                                                            sqlitedb = self.sqlitedb)

                
            ################################ Step 1. init
            self.script += "### Step 1: init db\n\n"
            self.script += """
{trino_cmd}init \\
\t--gene_trans_map {trans_map} \\
\t--transcript_fasta {trans_fa} \\
\t--transdecoder_pep  {pep_fa}\n\n""".format(trino_cmd = trino_cmd_sqlite, 
                                               trans_map = self.sample_data[sample]["gene_trans_map"], 
                                               trans_fa  = self.sample_data[sample]["transcripts.fasta.nucl"], 
                                               pep_fa    = self.sample_data[sample]["fasta.prot"])
            
            
            ################################ Step 2. Load
            self.script += "### Step 2: Load blast reports\n\n"
            self.script += """
{trino_cmd}LOAD_swissprot_blastp {blastp} \n
{trino_cmd}LOAD_swissprot_blastx {blastx} \n

""".format(trino_cmd = trino_cmd_sqlite,
           blastp  = self.sample_data[sample]["blast.prot"],
           blastx  = self.sample_data[sample]["blast.nucl"])

            if "hmmscan.prot" in self.sample_data[sample]:
                self.script += "{trino_cmd}LOAD_pfam {pfam}\n\n".format(trino_cmd=trino_cmd_sqlite,
                                                                        pfam=self.sample_data[sample]["hmmscan.prot"])

            if "rnammer" in self.sample_data[sample]:
                self.script += "{trino_cmd}LOAD_rnammer {rnammer}\n\n".format(trino_cmd=trino_cmd_sqlite,
                                                                              rnammer=self.sample_data[sample]["rnammer"])

            ################################ Step 4. Report
            self.script += "### Step 4: Create report\n\n"
            self.script += """
{trino_cmd}report \\
\t{redirects} > {dir}{file} \n\n""".format(trino_cmd = trino_cmd_sqlite,
                                 redirects = self.get_redir_parameters_script(),
                                 dir=use_dir, 
                                 file=output_basename)
            
            # Store results to fasta and assembly slots:
            self.sample_data[sample]["trino.rep"] = "%s%s" % (sample_dir, output_basename)
            
            self.stamp_file(self.sample_data[sample]["trino.rep"])

            ################################ Step 5. delete sqlitedb
            # If requested copy, create copy and change active sqlitedb
            if "cp_sqlitedb" in self.params and "rm_sqlitedb" in self.params:
                self.script += "rm -rf {dest}\n\n".format(dest=os.path.join(use_dir,os.path.basename(self.params["sqlitedb"])))
                        
            # Wrapping up function. Leave these lines at the end of every iteration:
            # Sees to copying local files to final destination (and other stuff)
            self.local_finish(use_dir,sample_dir)

            self.create_low_level_script()
