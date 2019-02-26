# -*- coding: UTF-8 -*-
""" 
``fasta_splitter``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for splitting `fasta` files into parts, using ``fasta-splitter.pl``.

Convenient for parallelizing processes on the cluster. You can take a project wide fasta file (such as a transcriptome), split it into sub-fasta files, and run various processes on the sub-files.

The parts can then be combined with ``merge_table`` module, which can concatenate any type of file.

The module ships with ``fasta-splitter.pl``  version 0.2.6, 2017-08-01

Usage::

    Usage: fasta-splitter [options] <file>...
        Options:
            --n-parts <N>        - Divide into <N> parts
            --part-size <N>      - Divide into parts of size <N>
            --measure (all|seq|count) - Specify whether all data, sequence length, or
                                   number of sequences is used for determining part
                                   sizes ('all' by default).
            --line-length        - Set output sequence line length, 0 for single line
                                   (default: 60).
            --eol (dos|mac|unix) - Choose end-of-line character ('unix' by default).
            --part-num-prefix T  - Put T before part number in file names (def.: .part-)
            --out-dir            - Specify output directory.
            --nopad              - Don't pad part numbers with 0.
            --version            - Show version.
            --help               - Show help.

You can't use the ``--part-size`` method, since it will end up in an unknown number of files, which is not defined in Neat-Seq Flow.

Please do **not** use the ``--nopad`` parameter. There is no reason to...

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A `fasta` file in one of the following slots (scope = "project"):

    * ``sample_data["project_data"]["fasta.nucl"]``
    * ``sample_data["project_data"]["fasta.prot"]``

* A `fasta` file in one of the following slots (scope = "sample"):

    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``



Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slots:
        
    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``

* For sample scope, the original sample list will be overridden with the new sample list.


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "nucl|prot", "The type of fasta file to split"
    "redirects --n-parts", "", "Number of fragments"
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    split_fasta1:
        module:         fasta_splitter
        base:           Trinity1
        script_path:    
        type:           nucl
        redirects:
            --n-parts:      4
            --measure:      seq

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`http://kirill-kryukov.com/study/tools/fasta-splitter/`_



"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"

class Step_fasta_splitter(Step):

    auto_redirs = "--out-dir --nopad --version --help".split(" ")

    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        if "type" not in self.params:
            raise AssertionExcept("Please supply a 'type' parameter: 'nucl' or 'prot'")
        if self.params["type"] not in ["nucl","prot"]:
            raise AssertionExcept("'type' parameter must be 'nucl' or 'prot'")
        if "scope" not in self.params:
            raise AssertionExcept("Please supply a 'scope' parameter: 'sample' or 'project'")
        if self.params["scope"] not in ["sample","project"]:
                raise AssertionExcept("'scope' parameter must be 'sample' or 'project'")

        # If user did not pass script_path, use the one packaged with the module
        if not self.params["script_path"]:
            cur_dir = os.path.dirname(os.path.realpath(__file__))
            self.params["script_path"] = "perl " + os.sep.join([cur_dir,"fasta-splitter.pl"])

        if "--n-parts" not in self.params["redir_params"]:
            raise  AssertionExcept("The module supports only the --n-parts version, since the number of chunks "
                                   "must be known at run-time! Please include the number of parts in the redirects")

        try:
            self.params["redir_params"]["--n-parts"] = int(self.params["redir_params"]["--n-parts"])
        except TypeError as type_err:
            raise AssertionExcept("--n-parts must be an integer, not a list or dictionary")
        except ValueError as val_err:
            raise AssertionExcept("--n-parts must be an integer or a string interpretable as an integer!")

        if "--measure" in self.params["redir_params"] and \
                self.params["redir_params"]["--measure"] not in ["all","seq","count"]:
            raise  AssertionExcept("--measure must be one of all, seq and count.")
        if "--part-num-prefix" in self.params["redir_params"]:
            if not isinstance(self.params["redir_params"]["--part-num-prefix"], str):
                raise  AssertionExcept("--part-num-prefix must be a string.")
        else:
            self.params["redir_params"]["--part-num-prefix"] = ".part-"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Check that scope exists and fasta.nucl or fasta.prot exist
        # Check that subsample_num exists
        # Check that type exists

        self.params["type"] = "fasta.{type}".format(type=self.params["type"])
        if self.params["scope"]=="project":
            if self.params["type"] not in self.sample_data["project_data"]:
                raise AssertionExcept("{type} does not exist in project".format(type=self.params["type"]))
        else:
            for sample in self.sample_data["samples"]:
                if self.params["type"] not in self.sample_data[sample]:
                    raise AssertionExcept("{type} does not exist in sample".format(type=self.params["type"]), sample)
                
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        if self.params["scope"] == "project":

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            self.script = ""
            # Get constant part of script:
            self.script += self.get_script_const()
            self.script += "--out-dir {dir} \\\n\t".format(dir=use_dir)
            self.script += "{fasta} \n\n".format(fasta=self.sample_data["project_data"][self.params["type"]])

            subsample_num = self.params["redir_params"]["--n-parts"]

            input_fn = os.path.basename(self.sample_data["project_data"][self.params["type"]])
            file_info = os.path.splitext(input_fn)
            num_of_digits = len(str(subsample_num))

            try:
                sample_list = ["subsample{num:0={res}}".format(num=num,res=num_of_digits)
                               for num
                               in range(1,int(subsample_num)+1)]
            except ValueError:
                raise AssertionExcept("'subsample_num' must be an integer")
            self.stash_sample_list(sample_list)

            # Creating data container for subsamples:
            for i in range(1, int(subsample_num)+1):             #self.sample_data["samples"]:
                # Formatting i to length of max i
                part_i = "{num:0={res}}".format(num=i,res=num_of_digits)
                # Formatting sample name (same as in loop above
                sample = "subsample{part}".format(part=part_i)
                self.sample_data[sample] = dict()
                self.sample_data[sample][self.params["type"]] = \
                    "{use_dir}{main}{prefix}{part}{ext}".format(use_dir=self.base_dir,
                                                                part=part_i,
                                                                main=file_info[0],
                                                                prefix=self.params["redir_params"]["--part-num-prefix"],
                                                                ext=file_info[1])
                # Stamping the files takes a long time. Cancelling for the time being
                # self.stamp_file(self.sample_data[sample][self.params["type"]])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                    
        else:  # self.params["scope"] == "sample"
        
            # raise AssertionExcept("Not defined yet...")
            # Each iteration must define the following class variables:
                # spec_script_name
                # script
            new_sample_list = list()
            new_sample_dict = dict()
            
            for sample in self.sample_data["samples"]:

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                self.script = ""
                # Get constant part of script:
                self.script += self.get_script_const()
                self.script += "--out-dir {dir} \\\n\t".format(dir=use_dir)
                self.script += "{fasta} \n\n".format(fasta=self.sample_data[sample][self.params["type"]])

                subsample_num = self.params["redir_params"]["--n-parts"]

                input_fn = os.path.basename(self.sample_data[sample][self.params["type"]])
                file_info = os.path.splitext(input_fn)
                num_of_digits = len(str(subsample_num))

                try:
                    subsample_list = ["{sample}.subsample{num:0={res}}".format(num=num, res=num_of_digits,sample=sample)
                                   for num
                                   in range(1, int(subsample_num) + 1)]
                except ValueError:
                    raise AssertionExcept("'subsample_num' must be an integer")
                new_sample_list.extend(subsample_list)

                # Creating data container for subsamples:
                for i in range(1, int(subsample_num) + 1):  # self.sample_data["samples"]:
                    # Formatting i to length of max i
                    part_i = "{num:0={res}}".format(num=i, res=num_of_digits)
                    # Formatting sample name (same as in loop above
                    subsample = "{sample}.subsample{part}".format(part=part_i, sample=sample)
                    self.sample_data[subsample] = dict()
                    self.sample_data[subsample][self.params["type"]] = \
                        "{use_dir}{main}{prefix}{part}{ext}".format(use_dir=sample_dir,
                                                                  part=part_i,
                                                                  prefix=self.params["redir_params"]["--part-num-prefix"],
                                                                  main=file_info[0],
                                                                  ext=file_info[1])
                    self.sample_data[subsample]["grouping"] = dict()
                    self.sample_data[subsample]["grouping"]["source"] = sample
                    self.sample_data[subsample]["type"] = self.determine_sample_types(subsample,self.sample_data[subsample])

                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir,sample_dir)
                self.create_low_level_script()

            self.sample_data["samples"] = new_sample_list
            self.stash_sample_list(new_sample_list)
