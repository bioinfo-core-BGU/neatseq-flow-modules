# -*- coding: UTF-8 -*-
""" 
``prokka_old``  :sup:`*`
-------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for running prokka:

Prokka is executed on the contigs stored in sample_data.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A nucleotide fasta file in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* If ``scope`` is set to ``sample``:

    * Puts output predicted protein sequences (faa file) in:
        ``sample_data[<sample>]["fasta.prot"]``

    * Puts output predicted protein genomic sequences (fna file) in:
        ``sample_data[<sample>]["fasta.nucl"]``

    * Puts the annotation file (gff) in:
        ``sample_data[<sample>]["gff"]``

    * Stores the prokks dir in:
        ``sample_data[<sample>]["prokka.dir"]``

* If ``scope`` is set to ``project``:

    * Puts output predicted protein sequences (faa file) in:
        ``sample_data["fasta.prot"]``

    * Puts output predicted protein genomic sequences (fna file) in:
        ``sample_data["fasta.nucl"]``

    * Puts the annotation file (gff) in:
        ``sample_data["gff"]``

    * Stores the prokks dir in:
        ``sample_data["prokka.dir"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "generate_GFF_dir", "empty", "Create a dir with links to the gff files for use downstream by others. Only relevant when ``scope=='sample'``"
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you set values to ``--locustag``, ``--genus``, ``--species`` and ``--strain``, these will hold for all the samples, and will be passed as-is to the scripts.

If you pass the parameters without setting their values, the values will be set to the sample names (or to the project name, when ``scope == 'project'``).


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    prokka1:
        module: prokka_old
        base: spades1
        script_path: /path/to/prokka
        qsub_params:
            -pe: shared 20
        generate_GFF_dir: 
        scope: sample
        redirects:
            --cpus: 20
            --fast: 
            --force:
            --genus: Staphylococcus
            --metagenome: 
            --strain: 

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seemann, T., 2014. **Prokka: rapid prokaryotic genome annotation**. *Bioinformatics*, 30(14), pp.2068-2069.
"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_prokka_old(Step):
    
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".prokka.out"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """


        sample_has_nucl = project_has_nucl = False
        if "scope" not in self.params:
            # If all samples have fasta.nucl:
            if all(map(lambda x: "fasta.nucl" in self.sample_data[x], self.sample_data["samples"])):
                sample_has_nucl = True
            if "fasta.nucl" in self.sample_data:
                project_has_nucl = True
            if sample_has_nucl and project_has_nucl:
                raise AssertionExcept("Both sample and project fasta exists. You must specify 'scope'")
            elif sample_has_nucl:
                self.params["scope"] = "sample"
            elif project_has_nucl:
                self.params["scope"] = "project"
            else:
                raise AssertionExcept("No fasta exists in either samples or project!")
            
            
        if self.params["scope"] == "sample":
            # Assert that all samples have nucleotide fasta files:
            for sample in self.sample_data["samples"]:    
                try:
                    self.sample_data[sample]["fasta.nucl"]
                except KeyError:
                    raise AssertionExcept("Sample does not have a fasta file\n", sample)
        elif self.params["scope"] == "project":
            try:
                self.sample_data["fasta.nucl"]
            except KeyError:
                raise AssertionExcept("Project does not have a fasta file\n")
        
        
        
            
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        if self.params["scope"] == "sample":
            if "generate_GFF_dir" in self.params.keys():
                
                gff_dir = "%sgff_file_links%s" % (self.base_dir, os.sep)
                
                self.script = "# Making dir with links to gff files\n\n"
                self.script += "mkdir %s\n\n" % gff_dir
                
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    link_name = "%s%s" % (gff_dir, os.path.basename(self.sample_data[sample]["gff"]))
                    self.script += "ln -s %s %s\n\n" % (self.sample_data[sample]["gff"], link_name)
        else:  # Scope = project
            pass
    
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
            
        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                
                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample])
                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                # Setting "--locustag","--genus","--species","--strain" to sample if they are passed but not set.
                for element in ["--locustag","--genus","--species","--strain"]:
                    if element in self.params["redir_params"].keys() and not self.params["redir_params"][element]:
                        self.params["redir_params"][element] = sample
                
                self.script += self.get_script_const()   
                self.script += "--outdir %s \\\n\t"  % sample_dir
                self.script += "--prefix %s \\\n\t"  % (sample + self.file_tag)

                self.script += "%s \n\n"           % (self.sample_data[sample]["fasta.nucl"])


                # Store results to fasta and assembly slots:
                self.sample_data[sample]["fasta.prot"]            = sample_dir + sample + self.file_tag + ".faa"
                self.sample_data[sample]["fasta.nucl"]            = sample_dir + sample + self.file_tag + ".fna"
                self.sample_data[sample]["gff"]            = sample_dir + sample + self.file_tag + ".gff"
                self.sample_data[sample]["prokka.dir"]    = sample_dir 



                
                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                            
                self.create_low_level_script()
                        
        else:  # scope = project
            
                            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)


            # Setting "--locustag","--genus","--species","--strain" to sample if they are passed but not set.
            for element in ["--locustag","--genus","--species","--strain"]:
                if element in self.params["redir_params"].keys() and not self.params["redir_params"][element]:
                    self.params["redir_params"][element] = self.sample_data["Title"]
            
            self.script += self.get_script_const()   
            self.script += "--outdir %s \\\n\t"  % self.base_dir
            self.script += "--prefix %s \\\n\t"  % (self.sample_data["Title"] + self.file_tag)

            self.script += "%s \n\n"           % (self.sample_data["fasta.nucl"])


            # Store results to fasta and assembly slots:
            self.sample_data["fasta.prot"]  = self.base_dir + self.sample_data["Title"] + self.file_tag + ".faa"
            self.sample_data["fasta.nucl"]  = self.base_dir + self.sample_data["Title"] + self.file_tag + ".fna"
            self.sample_data["gff"]         = self.base_dir + self.sample_data["Title"] + self.file_tag + ".gff"
            self.sample_data["prokka.dir"]  = self.base_dir 
                
                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            self.create_low_level_script()
                    
                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        pass