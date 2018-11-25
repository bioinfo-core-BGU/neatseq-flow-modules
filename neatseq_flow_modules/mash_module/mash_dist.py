# -*- coding: UTF-8 -*-
""" 
``mash_dist`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

    
Requires:
~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data["fasta.nucl"]``

* OR fastq files in one of the following slots (merge fastq files first with mash_sketch or otherwise):

    * ``sample_data[<sample>]["fastq"]``
    * ``sample_data["fastq"]``

* OR sketch files in one of the following slots:

    * ``sample_data[<sample>]["msh.fastq"]``
    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data["msh.fastq"]``
    * ``sample_data["msh.fasta"]``
    
    
Output:
~~~~~~~~~~~~~

* puts 'msh' output files in the following slots for (scope=sample):

    * ``sample_data[<sample>]["msh.fasta"]``
    * ``sample_data[<sample>]["msh.fastq"]``

* puts 'msh' output files in the following slots for (scope=project and scope=all_samples):

    * ``sample_data[<sample>]["mash.dist.table"]``
    * ``sample_data["mash.dist.table"]``
    

Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10
    
    "reference", "", "A block including 'path' or 'scope', 'type' and optionally 'msh'"
    "query", "", "A block including 'scope' (sample, project or all_samples), 'type' and optionally 'msh'"


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. External reference. Sample-wise fastq files.
    Returns table of mash dist of sample against external reference. One table per sample

::

    dist:
        module:         mash_dist
        base:           [sketch_proj,sketch_smp]
        script_path:    "{Vars.paths.mash} dist"
        reference:
            path:   /path/to/ref1
        query:
            scope:          sample
            type:           fastq
            msh:

2. Project mashed fasta reference. Sample mashed fastq query
    Returns table of mash dist of sample against project reference. One table per sample

::

    dist:
        module:         mash_dist
        base:           [sketch_proj,sketch_smp]
        script_path:    "{Vars.paths.mash} dist"
        reference:
            scope:      project
            type:       fasta
            msh:
        query:
            scope:      sample
            type:       fastq
            msh:

3. Project mashed reference. Project mashed fastq query
    Returns table of mash dist of project sketch against project sketch.
    One table for the whole project.

    If the project sketch is built from sample sketches, as is created by ``mash_sketch`` using ``scope=project`` and
    ``src_scope=sample``, the result will be an all-agianst-all mash dist table.

::

    dist:
        module:         mash_dist
        base:           [sketch_proj,sketch_smp]
        script_path:    "{Vars.paths.mash} dist"
        reference:
            scope:      project
            type:       fastq
            msh:
        query:
            scope:      project
            type:       fastq
            msh:

4. Project mashed fastq reference. Sample mashed fastq query
    Returns table of mash dist of project sketch against teach sample sketch. One table per sample.

::

    dist: 
        module:         mash_dist
        base:           [sketch_proj,sketch_smp]
        script_path:    "{Vars.paths.mash} dist"
        reference:
            scope:      project
            type:       fastq
            msh:
        query:
            scope:      sample
            type:       fastq
            msh:


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ondov, Brian D., et al. **Mash: fast genome and metagenome distance estimation using MinHash** *Genome biology*, 17.1 (2016): 132.

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step, AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_mash_dist(Step):
    
    auto_redirs = "-l".split(" ")

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".mash.tbl"

        if not (set(self.params) & {"reference" ,"query"}):
            raise AssertionExcept("You must pass 'reference' and 'query' blocks of params. See help")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # TODO: Check that files exist for the analysis
        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function
            
        """

        if self.params["query"]["scope"] in ["all_samples", "project"]:
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()

    def build_scripts_bysample(self):
        """ Script building function for sample-level BLAST
            
        """

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Setting reference sources:
            if "path" in self.params["reference"]:
                ref_path = self.params["reference"]["path"]
            else:
                if "type" in self.params["reference"] and self.params["reference"]["type"] == "fasta":
                    type2use = "fasta"
                else:
                    type2use = "fastq"

                if "scope" in self.params["reference"] and self.params["reference"]["scope"]=="project":
                    if "msh" in self.params["reference"]:
                        # if msh2use not in self.sample_data[sample]:
                        try:
                            ref_path = self.sample_data["project_data"]["msh." + type2use]
                        except KeyError:
                            raise AssertionExcept(
                                "No {type} mash sketch (in reference) for project".format(type=type2use))

                    else:
                        if type2use == "fasta":
                            type2use = "fasta.nucl"
                        try:
                            ref_path = self.sample_data["project_data"][type2use]
                        except KeyError:
                            raise AssertionExcept("No {type} file (in reference) for project".format(type=type2use))
                else:
                    if "msh" in self.params["reference"]:
                        # if msh2use not in self.sample_data[sample]:
                        try:
                            ref_path = self.sample_data[sample]["msh." + type2use]
                        except KeyError:
                            raise AssertionExcept("No {type} mash sketch (in reference) for sample".format(type=type2use), sample)

                    else:
                        if type2use == "fasta":
                            type2use = "fasta.nucl"
                        try:
                            ref_path = self.sample_data[sample][type2use]
                        except KeyError:
                            raise AssertionExcept("No {type} file (in reference) for sample".format(type=type2use), sample)

            # Setting query sources:
            if "type" in self.params["query"] and self.params["query"]["type"] == "fasta":
                type2use = "fasta"
            else:
                type2use = "fastq"
            if "msh" in self.params["query"]:

                try:
                    query_path = self.sample_data[sample]["msh." + type2use]
                except KeyError:
                    raise AssertionExcept("No {type} mash sketch (in query) for sample".format(type=type2use), sample)

            else:
                if type2use == "fasta":
                    type2use = "fasta.nucl"
                try:
                    query_path = self.sample_data[sample][type2use]
                except KeyError:
                    raise AssertionExcept("No {type} file (in query) for sample".format(type=type2use), sample)

            output_filename = sample + self.file_tag
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            self.script += self.get_script_const()
            self.script += "%s \\\n\t" % ref_path
            self.script += "%s \\\n\t" % query_path
            self.script += "> %s \n\n" % (use_dir + output_filename)

            # Store results table
            self.sample_data[sample]["mash.dist.table"] = (sample_dir + output_filename)
            self.stamp_file(self.sample_data[sample]["mash.dist.table"])

            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()
            
    def build_scripts_byproject(self):
        """ Script building function for project-level BLAST

        """

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        # Setting reference sources:
        if "path" in self.params["reference"]:
            ref_path = self.params["reference"]["path"]
        else:
            if "type" in self.params["reference"] and self.params["reference"]["type"] == "fasta":
                type2use = "fasta"
            else:
                type2use = "fastq"

            if "msh" in self.params["reference"]:
                # if msh2use not in self.sample_data[sample]:
                try:
                    ref_path = self.sample_data["project_data"]["msh." + type2use]
                except KeyError:
                    raise AssertionExcept("No {type} mash sketch (in reference) for project".format(type=type2use))

            else:
                if type2use == "fasta":
                    type2use = "fasta.nucl"
                try:
                    ref_path = self.sample_data["project_data"][type2use]
                except KeyError:
                    raise AssertionExcept("No {type} file (in reference) for project".format(type=type2use))

        # Setting query sources:
        if "type" in self.params["query"] and self.params["query"]["type"] == "fasta":
            type2use = "fasta"
        else:
            type2use = "fastq"
        if "msh" in self.params["query"]:

            if self.params["query"]["scope"] == "project":
                try:
                    query_path = self.sample_data["project_data"]["msh." + type2use]
                except KeyError:
                    raise AssertionExcept("No {type} mash sketch (in query) for project".format(type=type2use))
            else: # scope = "all_samples"
                try:
                    query_path = " ".join(
                        [self.sample_data[sample]["msh." + type2use]
                         for sample
                         in self.sample_data["samples"]])
                except KeyError:
                    raise AssertionExcept("A sample is missing a {type} mash sketch".format(type=type2use))
        else:
            if type2use == "fasta":
                type2use = "fasta.nucl"
            if self.params["query"]["scope"] == "project":
                try:
                    query_path = self.sample_data["project_data"][type2use]
                except KeyError:
                    raise AssertionExcept("No {type} file (in query) for project".format(type=type2use))
            else: # scope = "all_samples"
                try:
                    query_path = " ".join(
                        [self.sample_data[sample][type2use]
                         for sample
                         in self.sample_data["samples"]])
                except KeyError:
                    raise AssertionExcept("A sample is missing a {type} mash sketch".format(type=type2use))

        output_filename = self.sample_data["Title"] + self.file_tag

        self.script += self.get_script_const()
        self.script += "%s \\\n\t" % ref_path
        self.script += "%s \\\n\t" % query_path
        self.script += "> %s \n\n" % (use_dir + output_filename)

        # Store results table
        self.sample_data["project_data"]["mash.dist.table"] = (self.base_dir + output_filename)
        self.stamp_file(self.sample_data["project_data"]["mash.dist.table"])

        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)
        self.create_low_level_script()

