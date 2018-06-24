# -*- coding: UTF-8 -*-
""" 
``mash_sketch`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Build mash sketches from sequence files.

Works in three modes:

* ``scope=sample``
    Builds a separate sketch for each sample
* ``scope=project`` and ``src_scope=sample``
    Builds a project wide sketch from sample sequence files. This can be used with ``mash_dist`` module to perform
    all-against-all comparisons.
* ``scope=project``
    Builds a sketch from project sequence files.

    
Requires:
~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``

* or fastq files in the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
* For ``scope = project``, uses project-wide files.
    
Output:
~~~~~~~~~~~~~

* puts 'msh' output files in the following slots for (scope=sample):

    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data[<sample>]["msh.fastq"]``

* puts 'msh' output files in the following slots for (scope=project):

    * ``sample_data["msh.fasta"]``
    * ``sample_data["msh.fastq"]``



Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10
    
    "scope", "project|sample", "The scope for which to build the sketch."
    "src_scope", "project|sample", "The scope from which to take the sequence files. Default - same as ``scope``"
    "type", "nucl|prot", "Use fastq or fasta files. By default, uses any that exist."


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Create sketch for each sample based on fastq files

::

    sketch_smp:
        module:         mash_sketch
        base:           trim_gal
        script_path:    "{Vars.paths.mash} sketch"
        scope:          sample
        type:           fastq
        rm_merged:
        qsub_params:
            -pe:        shared 10
        redirects:
            -m:         2
            -p:         10

1. Create project sketch for all samples' fastq files

::

    sketch_proj:
        module:         mash_sketch
        base:           merge1
        script_path:    "{Vars.paths.mash} sketch"
        src_scope:      sample
        scope:          project
        type:           fastq
        rm_merged:
        redirects:
            -m:         2
            -p:         10


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ondov, Brian D., et al. **Mash: fast genome and metagenome distance estimation using MinHash** *Genome biology*, 17.1 (2016): 132.


"""
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_mash_sketch(Step):
    
    auto_redirs = "-o -l".split(" ")

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".msh"

        if self.params["scope"] not in ["sample","project"]:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        if "src_scope" in self.params:
            if self.params["src_scope"] not in ["sample", "project"]:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
            if self.params["src_scope"] == "project" and self.params["scope"] == "sample":
                raise AssertionExcept("Project 'src_scope' not defined for 'scope' sample.")

        else:
            self.params["src_scope"] = self.params["scope"]

        if "type" not in self.params:
            self.params["type"] = ["fastq","fasta"]
        else:
            if isinstance(self.params["type"], str):
                self.params["type"] = [self.params["type"]]

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # TODO: Check that files exist for the analysis
        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.script = ""
        if self.params["scope"]=="project" and self.params["src_scope"] == "sample":
            for type in list(set(self.params["type"]) & {"fastq", "fasta"}):

                # Create script only if there are files in files4mashing_lists
                if self.files4mashing_lists[type]:
                    # Define output filename 
                    output_filename = "".join([self.base_dir , self.sample_data["Title"] , ".", type])
                    # Note: The 'msh' tag in the output file get added automatically by mash

                    self.script += self.get_script_const()
                    self.script += "-o %s \\\n\t" % output_filename
                    self.script += "{msh_files} \n\n".format(msh_files=" \\\n\t".join(self.files4mashing_lists[type]))
                    
                    if "rm_merged" in self.params:  # and len(type_list)>1:
                        self.script += "# Removing temporary merged files\n"
                        for sample in self.files2remove:
                            for file in self.files2remove[sample]:
                                self.script += "rm -rf %s \n\n" % file

                    # Store msh file:
                    self.sample_data["msh." + type] = (output_filename + ".msh")
                    self.stamp_file(self.sample_data["msh." + type])

    def build_scripts(self):
        """ This is the actual script building function
            
        """

        if self.params["src_scope"] == "project":
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
        """

        self.files4mashing_lists = dict()
        self.files4mashing_lists["fasta"] = list()
        self.files4mashing_lists["fastq"] = list()
        self.files2remove = dict()  # Will contain list of temporary samples to remove at the end

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            self.files2remove[sample] = list()

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
            filetype_lists = {"fastq": list(set(self.sample_data[sample].keys()) & {"fastq.F", "fastq.R", "fastq.S"}),
                              "fasta": list(set(self.sample_data[sample].keys()) & {"fasta.nucl"})}

            # Loop over all types requested by user in 'type'
            for filetype in list(set(self.params["type"]) & {"fastq", "fasta"}):

                # If no files of type filetype exist, move on to next type
                if not filetype_lists[filetype]:
                    continue

                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,sample,filetype])
                self.script = ""

                # If there is only one sequence file (i.e. fastq.S), use directly
                if len(filetype_lists[filetype]) == 1:  # Only one file
                    input_filename = self.sample_data[sample][filetype_lists[filetype][0]]
                    final_input_filename = input_filename
                    # self.script += "# Single file in filetype %s. Not doing concatenation\n\n" % filetype
                # If there are more than one sequence file, concatenate.
                else:
                    input_filename = sample + ".merged." + filetype
                    # final_input_filename = sample_dir + sample + ".merged." + filetype
                    self.script += "cat \\\n\t"
                    for spec_type in filetype_lists[filetype]:
                        self.script += "%s \\\n\t" % self.sample_data[sample][spec_type]
                    self.script += "> %s\n\n" % (use_dir + input_filename)
                    self.files2remove[sample].append(use_dir + input_filename)
                    # This is done in here because sometimes the script is not finalized and 
                    # stamped files get carried from one script to the next.
                    self.stamp_file(sample_dir + input_filename)

                # Store the merged input to sketch:
                self.sample_data[sample][filetype] = sample_dir + input_filename

                # If scope is project, add file to files4mashing but don't do mashing.
                # Sketching itself is done in wrapping up after all merges are complete
                if self.params["scope"] == "project":
                    self.files4mashing_lists[filetype].append(sample_dir + input_filename)

                # If scope is sample, add code for actual mashing
                if self.params["scope"] == "sample":
                    # Define output filename 
                    output_filename = "%s.%s" % (sample, filetype)  # The 'msh' tag get added automatically by mash

                    self.script += self.get_script_const()
                    self.script += "-o %s \\\n\t" % (use_dir + output_filename)
                    self.script += "%s \n\n" % input_filename
                    
                    if "rm_merged" in self.params:
                        for file in self.files2remove[sample]:
                            self.script += "rm -rf %s \n\n" % file

                    # Store msh  file:
                    self.sample_data[sample]["msh." + filetype] = (sample_dir + output_filename + ".msh")
                    self.stamp_file(self.sample_data[sample]["msh." + filetype])

                # If self.script was created, create the script file
                # e.g. if scope is project, src_scope is sample and there is only one sequence file,
                # no script will be created!
                if self.script:
                    # Wrapping up function. Leave these lines at the end of every iteration:
                    self.local_finish(use_dir,sample_dir)
                    self.create_low_level_script()
            
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        type_counter = 0 # Counts the number of types for which info was found. Fails if no files are found

        for filetype in list(set(self.params["type"]) & {"fastq","fasta"}):
            # Name of specific script:
            self.spec_script_name = "_".join([self.step, self.name, self.sample_data["Title"], filetype])
            self.script = ""

            filetype_lists = {"fastq": list(set(self.sample_data.keys()) & {"fastq.F", "fastq.R", "fastq.S"}),
                              "fasta": list(set(self.sample_data.keys()) & {"fasta.nucl"})}

            print filetype
            print filetype_lists
            # if filetype=="fastq":
            #     type_list = list(set(self.sample_data.keys()) & set(["fastq.F","fastq.R","fastq.S"]))
            # else:
            #     type_list = list(set(self.sample_data.keys()) & set(["fasta.nucl"]))
            # import pdb; pdb.set_trace()

            # Move on if no files exist for type
            if not filetype_lists[filetype]:
                continue
            # Count the number of types for which info was found. Fails if no files are found
            type_counter += 1

            # If there is only one file, use it
            # final_input_filename points to source file
            if len(filetype_lists[filetype]) == 1:  # Only one file
                input_filename = self.sample_data[filetype_lists[filetype][0]]
                final_input_filename = input_filename
            # If there is are more than one file, concatenate them
            # final_input_filename points to the new concatenated file
            else:
                input_filename       = use_dir       + self.sample_data["Title"] + ".merged." + filetype
                final_input_filename = self.base_dir + self.sample_data["Title"] + ".merged." + filetype
                self.script += "cat \\\n\t"
                for spec_type in filetype_lists[filetype]:
                    self.script += "%s \\\n\t" % self.sample_data[spec_type]
                self.script += "> %s\n\n" % input_filename


            # Define output filename 
            output_filename = "%s.%s" % (self.sample_data["Title"] , filetype)
            # The 'msh' tag get added automatically by mash

            self.script += self.get_script_const()
            self.script += "-o %s \\\n\t" % (use_dir + output_filename)
            self.script += "%s \n\n" % input_filename
            
            if "rm_merged" in self.params and len(filetype_lists[filetype])>1:
                self.script += "rm -rf %s \n\n" % input_filename

            # Store msh file:
            self.sample_data["msh." + filetype] = (self.base_dir + output_filename + ".msh")
            self.sample_data[filetype] = final_input_filename
            self.stamp_file(self.sample_data["msh." + filetype])
            self.stamp_file(self.sample_data[filetype])
            
            # Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
        
        if type_counter == 0:
            raise AssertionExcept("""\
No source data was found for type {types}. Check definitions of 'scope' and 'src_scope'""".
                                  format(types=", ".join(list(set(self.params["type"]) & {"fastq", "fasta"}))))
