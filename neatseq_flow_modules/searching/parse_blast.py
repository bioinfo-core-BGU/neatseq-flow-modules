# -*- coding: UTF-8 -*-
""" 
``parse_blast``
-----------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for running ``parse_blast_general.R``:

The ``parse_blast_general`` script will be available on github in the near future.

The program performs the following tasks:

1. It adds annotation to raw tabular BLAST output files,
2. filters the BLAST results by several possible fields,
3. selects the best hit for a group when passed a grouping field and
4. extracts the sequences equivalent to the alignments.



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Tabular BLAST result files in the following slots:

    * ``sample_data[<sample>]["blast.nucl|blast.prot"]`` if ``scope = sample``
    * ``sample_data["blast.nucl|blast.prot"]``           if ``scope = project``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the parsed report in:

    * ``self.sample_data[<sample>]["blast.parsed"]``  if ``scope = sample``
    * ``self.sample_data["blast.parsed"]``            if ``scope = project``

    

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "fasta2use", "``nucl|prot``", "If both nucl and prot BLAST reports exist, you have to specify which one to use with this parameter."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    parse_blast_table:
        module: parse_blast
        base: blst_table
        script_path: {Vars.paths.parse_blast}
        blast_merge_path: '{Vars.paths.compare_blast_parsed_reports} --variable evalue --full_txt_output'
        scope: sample
        redirects:
            --columns2keep: '"group name accession qseqid sallseqid evalue bitscore
                score pident coverage align_len"'
            --dbtable: {Vars.databases.gene_list.table}
            --group_dif_name: # See parse_blast_general.R documentation for how this is to be specified
            --max_evalue: 1e-7
            --merge_blast: qseqid
            --merge_metadata: # See parse_blast_general.R documentation for how this is to be specified
            --min_align_len: 30
            --min_coverage: 60
            --names: '"qseqid sallseqid qlen slen qstart qend sstart send length evalue
                bitscore score pident qframe"'
            --num_hits: 1


"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_parse_blast(Step):
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

        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                if not "blast" in self.sample_data:
                    raise AssertionExcept("There are no project BLAST results.\n")
            elif self.params["scope"]=="sample":
                # Checking all samples have a 'blast' file-type in sample_data
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    if not "blast" in self.sample_data[sample]:
                        raise AssertionExcept("There are no BLAST results.\n" , sample)
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
                
            
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.make_sample_file_index()   # see definition below
        
                
        try:
            self.params["blast_merge_path"]
        except KeyError:
            self.write_warning("You did not supply a blast_merge_path parameter. Will not merge blast reports...\n")
            self.script = ""
        else:
            self.script = "%s \\\n\t--blastind %s \\\n\t--output %s\n\n" % (self.params["blast_merge_path"], \
                                                                            self.sample_data["BLAST_files_index"], \
                                                                            self.base_dir + "merged_parsed_blast.tsv")
    
    def build_scripts(self):
        
        if self.params["scope"]=="project":
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
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
                if "fasta2use" in self.params.keys() and self.params["fasta2use"] in ("nucl","prot"):
                    fasta2use = self.params["fasta2use"]
                    # self.script += "--blast %s \\\n\t" % self.sample_data[sample]["blast"][fasta2use]
                else:
                    raise AssertionExcept("Project has both 'nucl' and 'prot' blast results. Select one by specifying the 'fasta2use' parameter.")
            elif "blast.nucl" in self.sample_data:
                fasta2use = "nucl"
            elif "blast.prot" in self.sample_data:
                fasta2use = "prot"
            else:
                raise AssertionExcept("No BLAST Results defined\n")

            self.script += self.get_script_const()
            self.script += "--blast %s \\\n\t" % self.sample_data["blast." + fasta2use]
            
            # FASTA Extraction
            if "extract_fasta" in self.params:
                try:
                    self.script += "--fasta2extract %s \\\n\t" % self.sample_data["fasta." + fasta2use]
                except keyError:
                    raise AssertionExcept("In order to extract the fasta sequences, you need to have a project wide fasta file defined with the same type as the blast type.")


            self.script += "--output %s\n\n" % os.sep.join([use_dir,output_filename])
            
            
            # Store BLAST result file:
            self.sample_data["blast.parsed"] = "".join([self.base_dir, output_filename])
            self.sample_data["blast.parsed." + fasta2use] = self.sample_data["blast.parsed"]
            self.stamp_file(self.sample_data["blast.parsed"])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                      
            
            self.create_low_level_script()
                        
                
        else:  # self.params["scope"]=="sample":
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

                # Define query and db files:
                # If db is defined by user, set the query to the correct 'fasta2use'
                # If both nucl and prot appear in blast results
                if "blast.nucl" in self.sample_data[sample] and "blast.prot" in self.sample_data[sample]:
                    if "fasta2use" in self.params.keys() and self.params["fasta2use"] in ("nucl","prot"):
                        fasta2use = self.params["fasta2use"]
                        # self.script += "--blast %s \\\n\t" % self.sample_data[sample]["blast"][fasta2use]
                    else:
                        raise AssertionExcept("Sample has both 'nucl' and 'prot' blast results. Select one by specifying the 'fasta2use' parameter." , sample)
                elif "blast.nucl" in self.sample_data[sample]:
                    fasta2use = "nucl"
                elif "blast.prot" in self.sample_data[sample]:
                    fasta2use = "prot"
                else:
                    raise AssertionExcept("No BLAST Results defined\n")

                # Define the actual script:    
                self.script += self.get_script_const()
                self.script += "--blast %s \\\n\t" % self.sample_data[sample]["blast." + fasta2use]
                
                # FASTA Extraction
                if "extract_fasta" in self.params:
                    try:
                        self.script += "--fasta2extract %s \\\n\t" % self.sample_data[sample]["fasta." + fasta2use]
                    except keyError:
                        raise AssertionExcept("In order to extract the fasta sequences, you need to have a fasta file defined with the same type as the blast type.", sample)

                self.script += "--output %s\n\n" % output_filename
                
                
                # Store BLAST result file:
                self.sample_data[sample]["blast.parsed"] = "".join([sample_dir , sample , self.file_tag])
                self.sample_data[sample]["blast.parsed." + fasta2use] = self.sample_data[sample]["blast.parsed"]
                self.stamp_file(self.sample_data[sample]["blast.parsed"])

                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                          
                
                self.create_low_level_script()
                        
                
   

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "parsed_BLAST_files_index.txt", "w") as index_fh:
            index_fh.write("#Sample\tBLAST_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["blast.parsed"]))
                
        self.sample_data["BLAST_files_index"] = self.base_dir + "parsed_BLAST_files_index.txt"
        
  