# -*- coding: UTF-8 -*-
"""
``bamtools``
-----------------------



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BAM file in one of the following slots:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data["project_data"]["bam"]``



Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Depending on the parameters:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["fastq.I"]``  (Interleaved fastq file)
    * ``sample_data[<sample>]["bed|fasta.nucl|json|pileup|sam|yaml"]``
    * ``sample_data[<sample>]["bai"]``
    * ``sample_data[<sample>]["bam.count"]``
    * ``sample_data[<sample>]["bam.coverage"]``
    * ``sample_data[<sample>]["bam.stats"]``
    * ``sample_data[<sample>]["bam.header"]``

For project scope, will go in project_data.

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "``project|sample", "The origin of the BAM file to work on."


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    bamtools_bwt1:
        module:             bamtools
        base:               sam_bwt1
        script_path:        {Vars.paths.bamtools}
        filter:             -mapQuality ">50"
        sort:
        random:
        count:
        index:
        coverage:

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""
import os
import sys
import re

from neatseq_flow.PLC_step import Step, AssertionExcept

__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_bamtools(Step):
    auto_redirs = "--in --out".split(" ")


    def step_specific_init(self):
        self.shell = "bash"  # Can be set to "bash" by inheriting instances

        self.bam_output = "filter random resolve revert sort".split(" ")
        self.non_bam_output = ["count coverage index stats header".split(" ")]
        self.not_supported = ["split"]
        self.specially_supported = ["merge convert"]

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if "scope" not in self.params:
            self.params["scope"] = "sample"
            self.write_warning("No scope defined. Using sample scope by default.", admonition="attention".upper())

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash
            if "bam" not in self.sample_data[sample]:
                raise AssertionExcept("No bam file found...", sample)

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
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)

            # outgtf = ".joined".join(os.path.splitext(os.path.basename(self.sample_data[sample]["gtf"])))
            active_file = self.sample_data[sample]["bam"]
            self.script = ""

            # If there are bam-creating commands:
            if list(set(self.params.keys()) & set(self.bam_output)):
                # Suffix is concatenation of bam-creating commands
                suffix = "." + ".".join(list(set(self.params.keys()) & set(self.bam_output)))
                outfile = suffix.join(os.path.splitext(os.path.basename(active_file)))

                first_cmd = True
                for cmd in self.params:
                    if cmd in self.bam_output:
                        if not first_cmd:
                            self.script += " \\\n| "
                        self.script += " \\\n\t".join([x for x in [self.params["script_path"],
                                                             cmd,
                                                             "-in " + active_file if first_cmd else None,
                                                             self.params[cmd]] if x is not None])

                        first_cmd = False
                self.script += " \\\n\t-out {bam}\n\n".format(bam=use_dir+outfile)

                self.sample_data[sample]["bam"] = "{dir}{bam}".format(dir=sample_dir,bam=outfile)
                self.stamp_file(self.sample_data[sample]["bam"])

                active_file = use_dir + outfile
            else:  # No bam creating steps. Creating local link to original bam file
                self.script += """\
##########
# Making local link to original bam file: (-f to force)
#----------
cp -fs {active_file} {here}

""".format(active_file=active_file,
           here=use_dir)
                active_file = use_dir + os.path.basename(active_file)


            for cmd in self.params:
                if cmd == "convert":
                    try:
                        re1 = re.search(pattern="-format\s+(\w+)",
                                        string=self.params[cmd])
                        totype = re1.group(1)
                    except AttributeError:
                        raise AssertionExcept("You must supply a -format flag in 'convert'")
                    if totype not in "bed fasta fastq json pileup sam yaml".split(" "):
                        raise AssertionExcept("-format must be one of: 'bed fasta fastq json pileup sam yaml'")
                    outfile = ".".join([os.path.basename(active_file),totype])

                    self.script += "\n\n"+" \\\n\t".join([x for x in [self.params["script_path"],
                                                                 cmd,
                                                                 "-in " + active_file,
                                                                 self.params[cmd],
                                                                 "-out " + use_dir + outfile] if x is not None])
                    if totype=="fasta":
                        totype = "fasta.nucl"
                    if totype=="fastq":
                        totype = "fastq.I"  # Interleaved fastq!
                    self.sample_data[sample][totype] = sample_dir + outfile
                    self.stamp_file(self.sample_data[sample][totype])

                # if cmd in self.non_bam_output:
                if cmd in ["count", "header", "stats", "coverage"]:
                    self.script += "\n\n" + " \\\n\t".join([x for x in [self.params["script_path"],
                                                                   cmd,
                                                                   "-in " + active_file,
                                                                   self.params[cmd],
                                                                   " > " + ".".join([active_file,cmd,"txt"])] if x is not None])
                    self.sample_data[sample]["bam."+cmd] = sample_dir + ".".join([active_file,cmd,"txt"])
                    self.stamp_file(self.sample_data[sample]["bam."+cmd])
                elif cmd in ["index"]:
                    self.script += "\n\n" + " \\\n\t".join([x for x in [self.params["script_path"],
                                                                   cmd,
                                                                   "-in " + active_file,
                                                                   self.params[cmd]] if x is not None])
                    self.sample_data[sample]["bai"] = sample_dir + ".".join([active_file,cmd,"txt"])
                    self.stamp_file(self.sample_data[sample]["bam."+cmd])
                else:
                    pass




            self.local_finish(use_dir, sample_dir)
            self.create_low_level_script()


