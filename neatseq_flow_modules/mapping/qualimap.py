# -*- coding: UTF-8 -*-
""" 
``qualimap`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running qualimap bamqc and multi-bamqc.

In the future it will include the other qualimap functionality

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* bam files in one of the following slots:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data["project_data"]["bam"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
.. * puts fastqc output files in the following slots:
..
..     * ``sample_data[<sample>]["fastqc_fastq.F_html"]``
..     * ``sample_data[<sample>]["fastqc_fastq.R_html"]``
..     * ``sample_data[<sample>]["fastqc_fastq.S_html"]``
..
.. * puts fastqc zip files in the following slots:
..
..     * ``sample_data[<sample>]["fastqc_fastq.F_zip"]``
..     * ``sample_data[<sample>]["fastqc_fastq.R_zip"]``
..     * ``sample_data[<sample>]["fastqc_fastq.S_zip"]``
            
 

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


::

    Qualimap_bamqc:
        module:         qualimap
        base:           Samtools
        script_path:    {Vars.paths.qualimap}
        scope:          sample
        mode:           bamqc

    Qualimap_multibamqc:
        module:         qualimap
        base:           Qualimap_bamqc
        script_path:    {Vars.paths.qualimap}
        scope:          project
        mode:           multi-bamqc



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_qualimap(Step):

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "fasqc"

        self.auto_redirs = "-d --data -outdir -outfile -bam".split(" ")

        if "mode" not in self.params:
            raise AssertionExcept("Please supply a qualimap mode")

        if self.params["mode"] not in ["bamqc", "multi-bamqc", "rnaseq"]:
            raise AssertionExcept("mode must be 'bamqc', 'rnaseq' or 'multi-bamqc'")

        if "scope" not in self.params:
            if self.params["mode"] == "multi-bamqc":
                self.params["scope"] = "project"
            else:
                self.params["scope"] = "sample"
            self.write_warning("No scope specified. Setting to '{scope}'".format(scope=self.params["scope"]))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            pass

        if self.params["mode"] in ["bamqc","rnaseq"]:
            for sample in sample_list:
                if "bam" not in self.sample_data[sample]:
                    raise AssertionExcept("No bam file defined for sample\n", sample)
        elif self.params["mode"] == "multi-bamqc":
            if "qualimap_files_index" not in self.sample_data["project_data"]:
                raise AssertionExcept("if mode is multi-bamqc, you have to have a 'qualimap bamqc' instance first!")
        else:
            pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        if self.params["scope"] == "sample":
            self.make_sample_file_index()   # see definition below

    def make_sample_file_index(self):
        """ Make file containing samples and target file names.
            This can be used by scripts called by create_spec_wrapping_up_script() to summarize the BLAST outputs.
        """

        with open(self.base_dir + "qualimap_files_index.txt", "w") as index_fh:
            # index_fh.write("Sample\tBLAST_report\n")
            for sample in self.sample_data["samples"]:  # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample, self.sample_data[sample]["qualimap"]))

        self.sample_data["project_data"]["qualimap_files_index"] = self.base_dir + "qualimap_files_index.txt"

    def build_scripts(self):
        """ This is the actual script building function
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:      # Getting list of samples out of samples_hash

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            annot2use = None
            for feat in ["gff","gtf","bed"]:
                if feat in self.params:
                    if self.params[feat] is None:
                        if feat in self.sample_data[sample]:
                            annot2use = self.sample_data[sample][feat]
                        elif feat in self.sample_data["project_data"]:
                            annot2use = self.sample_data["project_data"][feat]
                        else:
                            raise AssertionExcept("You passed an empty '{feat}' but none exists in project!".format(feat=feat))
                    else:
                        annot2use = self.params[feat]


            if self.params["mode"] == "bamqc":
                self.script = """
{setenv}
{const} bamqc \\
	{redir}
	-bam {bam} \\{feat}
	-outdir {outd} \\
	-outfile {outf}
            """.format(const=self.params["script_path"],
                       setenv=self.get_setenv_part(),
                       redir=self.get_redir_parameters_script().rstrip(),
                       bam=self.sample_data[sample]["bam"],
                       feat="" if not annot2use else "\n\t-gff " + annot2use + " \\",
                       outd=use_dir,
                       outf="{smp}_qualimap.report".format(smp=sample))

            elif self.params["mode"] == "multi-bamqc":
                self.script = """

{setenv}
{const} multi-bamqc \\
	{redir}
	--data {report} \\
	-outdir {outd} \\
	-outfile {outf}
                """.format(const=self.params["script_path"],
                           setenv=self.get_setenv_part(),
                           redir=self.get_redir_parameters_script().rstrip(),
                           report=self.sample_data[sample]["qualimap_files_index"],
                           outd=use_dir,
                           outf="{proj}_qualimap.report".format(proj=self.sample_data["Title"]))

            elif self.params["mode"] == "rnaseq":
                self.script = """

{setenv}
{const} rnaseq \\
	{redir}
	-bam {bam} \\{feat}
	-outdir {outd} \\
	-outfile {outf}
""".format(const=self.params["script_path"],
           setenv=self.get_setenv_part(),
           redir=self.get_redir_parameters_script().rstrip(),
           bam=self.sample_data[sample]["bam"],
           feat="" if not annot2use else "\n\t-gtf " + annot2use + " \\",
           outd=use_dir,
           outf="{proj}_qualimap.report".format(proj=self.sample_data["Title"]))

            self.sample_data[sample]["qualimap"] = sample_dir
            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()


#     def build_scripts_byproject(self):
#         # Each iteration must define the following class variables:
#         # spec_script_name
#         # script
#
#         # Name of specific script:
#         self.spec_script_name = self.set_spec_script_name()
#         self.script = ""
#
#         # This line should be left before every new script. It sees to local issues.
#         # Use the dir it returns as the base_dir for this step.
#         use_dir = self.local_start(self.base_dir)
#
#         if self.params["mode"] == "bamqc":
#             self.script = """
#
# {setenv}
# {const} {mode} \\
# 	{redir}
# 	-bam {bam} \\
# 	-outdir {outd} \\
# 	-outfile {outf}
#             """.format(const=self.params["script_path"],
#                        mode=self.params["mode"],
#                        setenv=self.get_setenv_part(),
#                        redir=self.get_redir_parameters_script().rstrip(),
#                        bam=self.sample_data["project_data"]["bam"],
#                        outd=use_dir,
#                        outf="{proj}_qualimap.report".format(proj=self.sample_data["Title"]))
#
#             self.sample_data["project_data"]["qualimap"] = self.base_dir
#
#
#         elif self.params["mode"] == "multi-bamqc":
#             self.script = """
#
# {setenv}
# {const} multi-bamqc \\
# 	{redir}
# 	--data {report} \\
# 	-outdir {outd} \\
# 	-outfile {outf}
#             """.format(const=self.params["script_path"],
#                        setenv=self.get_setenv_part(),
#                        redir=self.get_redir_parameters_script().rstrip(),
#                        report=self.sample_data["project_data"]["qualimap_files_index"],
#                        outd=use_dir,
#                        outf="{proj}_qualimap.report".format(proj=self.sample_data["Title"]))
#
#             self.sample_data["project_data"]["qualimap"] = self.base_dir
#
#         self.local_finish(use_dir, self.base_dir)
#         self.create_low_level_script()

