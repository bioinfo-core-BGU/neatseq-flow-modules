# -*- coding: UTF-8 -*-
"""
``samtools_new`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for executing samtools on a SAM or BAM file.

.. attention:: The module was tested on samtools 1.3

The samtools programs included in the module are the following:

* ``view`` to convert the SAM file to a BAM file
* ``sort`` to sort the BAM file
* ``index`` creates an index for the BAM file
* ``flagstat`` Runs flagstat on the BAM file
* ``stats`` Runs stats on the BAM file
* ``idxstats`` Runs idxstats on the BAM file
* ``fastq/a`` Converts a BAM or CRAM into either FASTQ or FASTA format depending on the command invoked.
* ``merge`` Merges sample bam files into single project bam file.

.. Note:: Order of samtools subprogram execution:

    The ``samtools`` programs are executed in the order above. It is up to you to have a sensible combination...


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* A SAM file in the following location:

    * ``sample_data[<sample>]["sam"]`` (for ``scope=sample``)
    * ``sample_data["sam"]`` (for ``scope=project``)

* Or a BAM file in:

    * ``sample_data[<sample>]["bam"]`` (for ``scope=sample``)
    * ``sample_data["bam"]`` (for ``scope=project``)

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Depending on the parameters, will put files in the following locations:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["bai"]``
    * ``sample_data[<sample>]["unfiltered_bam"]``
    * ``sample_data[<sample>]["unsorted_bam"]``
    * ``sample_data[<sample>]["bam.flagstat"]``
    * ``sample_data[<sample>]["bam.stats"]``
    * ``sample_data[<sample>]["bam.idxstats"]``

* If ``fastq`` was called, will also create the following files:

    * ``self.sample_data[<sample>]["fastq.F"]``
    * ``self.sample_data[<sample>]["fastq.R"]``
    * ``self.sample_data[<sample>]["fastq.S"]``

* If ``fasta`` was called, will also create the following files:

    * ``self.sample_data[<sample>]["fasta.F"]``
    * ``self.sample_data[<sample>]["fasta.R"]``
    * ``self.sample_data[<sample>]["fasta.S"]``


.. Note:: If ``scope`` is set to ``project``, the above mentioned output files will be created in the project
    scope, e.g. ``sample_data["project_data"]["stats"]``..

.. Note:: If ``merge`` is included, ``scope`` must be ``sample`` and the merged *bam* is located in ``sample_data["project_data"]["stats"]``..

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "project", "sample|project", "Scope of SAM/BAM top operate on. Defaults to ``sample``."
    "view", "*e.g.*: -buh  -q 30", "``samtools view`` parameters."
    "sort", "*e.g.*: -@ 20", "``samtools sort`` parameters."
    "index", "", "``samtools index`` parameters."
    "flagstat", "", "Leave empty. flagstat takes no parameters"
    "stats", "``samtools stats`` parameters", "Adds code for ``samtools stats``"
    "idxstats", "", "Adds code for ``samtools idxstats``"
    "fastq/a", "``samtools fastq/a`` parameters", "Adds code for ``samtools fastq/a``"
    "merge", "``*e.g.*: -R region``", "Adds code for ``samtools merge``, using the parameters supplied"
    "region", "", "A region to limit the ``view`` script to."
    "filter_by_tag", "*e.g.*: NM:i:[01]", "Filter BAM by one of the tags. Use an awk-compliant regular expression. In this example, keep only lines where the edit distance is 0 or 1. This is an experimental feature and should be used with caution..."
    "del_sam", "", "Remove SAM file"
    "del_unsorted", "", "Remove unsorted bam file."
    "type2use","sam|bam","Type of file to use. Must exist in scope"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sam_bwt1:
        module: samtools
        base: bwt1
        script_path: /path/to/samtools/bin/samtools
        qsub_params:
            -pe: shared 20
        view: -buh  -q 30 -@ 20 -F 4
        sort: -@ 20
        flagstat:
        index:
        stats: --remove-dups
        del_sam:
        del_unsorted:


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G. and Durbin, R., 2009. **The sequence alignment/map format and SAMtools**. *Bioinformatics*, 25(16), pp.2078-2079.

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept
import yaml
from pprint import pprint as pp


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_samtools_new(Step):
       

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        fastaq_actions = "dict faidx fqidx".split(" ")
        fastaq_actions = list(set(self.params) & set(fastaq_actions))
        if len(fastaq_actions) > 0:
            raise AssertionExcept("These tools are not supported by this module: {tools}".format(tools=", ".join(fastaq_actions)))

        if "type2use" in self.params:
            if not isinstance(self.params["type2use"],str) or self.params["type2use"] not in ["sam","bam"]:
                raise AssertionExcept("'type2use' must be either 'sam' or 'bam'")

        # Setting default scope to sample
        if "scope" not in self.params:
            self.params["scope"] = "sample"

        # for prog in "view sort index flagstat stats idxstats fastq fasta merge".split(" "):
        #     if prog in self.params and self.params[prog] is None:
        #         self.params[prog] = ""

        # Load YAML of file type stored in merge_file_types.yml
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "samtools_params.yml"), "r") as fileh:
            try:
                self.samtools_params = yaml.load("".join(fileh.readlines()), Loader=yaml.SafeLoader)
            except yaml.YAMLError as exc:
                if hasattr(exc, 'problem_mark'):
                    mark = exc.problem_mark
                    print("Error position: (%s:%s)" % (mark.line + 1, mark.column + 1))
                    print(mark.get_snippet())
                raise AssertionExcept("Error loading samtools params index 'samtools_params.yml'")
            except:
                raise AssertionExcept("Unknown error loading samtools params index 'samtools_params.yml")

        import json

        if "keep_output" not in self.params:
            self.params["keep_output"] = list(set(self.params.keys()) & set(self.samtools_params.keys()))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Set list of samples to go over. Either self.sample_data["samples"] for sample scope
        # or ["project_data"] for project scope
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
            if "merge" in list(self.params.keys()):
                raise AssertionExcept("project scope not defined for samtools merge")
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Check that a sam or bam exists
            sambam_actions = "addreplacerg calmd depad fasta fastq fixmate markdup merge phase index mpileup " \
                             "reheader bedcov collate depth flagstat idxstats quickcheck sort split stats view".split(" ")
            sambam_actions = list(set(self.params) & set(sambam_actions))
            if len(sambam_actions) > 0:
                if "bam" in self.sample_data[sample] and "sam" in self.sample_data[sample]:
                    if "type2use" in self.params:
                        self.file2use = self.params["type2use"]
                    else:
                        raise AssertionExcept(
                            "Both BAM and SAM file types exist. Specify which one to use with 'type2use'.\n",
                            sample)
                elif "bam" in self.sample_data[sample]:
                    self.file2use = "bam"
                elif "sam" in self.sample_data[sample]:
                    self.file2use = "sam"
                else:
                    raise AssertionExcept("Neither BAM nor SAM file exist for actions {actions}.\n".format(actions=", ".join(sambam_actions)), sample)


        if "type2use" in self.params:
            if self.params["type2use"] not in self.sample_data[sample]:
                raise AssertionExcept("No file of type '{type}' exists.".format(type=self.params["type2use"]),
                                      sample)


    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
#         # -------------- samtools merge ------------------
#         if "merge" in list(self.params.keys()):
#             sample_dir = self.make_folder_for_sample()
#
#             # Name of specific script:
#             self.spec_script_name = self.set_spec_script_name()
#             self.script = ""
#
#             # This line should be left before every new script. It sees to local issues.
#             # Use the dir it returns as the base_dir for this step.
#             use_dir = self.local_start(sample_dir)
#
#             outfile = self.sample_data["Title"] + ".merged.bam"
#
#             self.script += """\
# ###########
# # Running samtools merge
# #----------------
#
# {env_path} merge \\{params}
# \t{outfile} \\
# \t{infiles}
#
# """.format(env_path=self.get_script_env_path(),
#            infiles=" \\\n\t".join([self.sample_data[sample]["bam"] for sample in self.sample_data["samples"]]),
#            params="" if not self.params["merge"] else "\n\t" + self.params["merge"] + " \\",
#            outfile=use_dir + outfile)
#
#             self.sample_data["project_data"]["bam"] = sample_dir + outfile
#             self.stamp_file(self.sample_data["project_data"]["bam"])

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

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Will create all files in temp, and move only final version of each type to final location
            # Then, will remove temp
            temp_use_dir = use_dir+"temp"+os.sep
            try:
                os.makedirs(temp_use_dir)
            except FileExistsError:
                pass

            # active_file = self.sample_data[sample][self.file2use]
            active_files = dict(zip((self.file2use,),(self.sample_data[sample][self.file2use],)))
            active_type = self.file2use

            files2keep = list()

            # Starting off with local link to active file:
            self.script += """\
##########
# Making local link to original bam file: (-f to force)
#----------
cp -fs \\
\t{active_file} \\
\t{here}

""".format(active_file=active_files[active_type],
           here=temp_use_dir)

            active_files[active_type] = temp_use_dir + os.path.basename(active_files[active_type])
            self.sample_data[sample][active_type] = sample_dir + os.path.basename(active_files[active_type])

            for action in self.params:
                # This is to enable passing each tool more than once:
                # Repeated actions have digits appended to them
                action_numbered = action
                action = re.sub("\d+$","",action)

                if action not in self.samtools_params:
                    continue

                # Get action redirects (if exist)
                redirects = ""
                region = ""
                if isinstance(self.params[action_numbered],dict):
                    if "region" in self.params[action_numbered]:
                        region = self.params[action_numbered]["region"]
                    if "redirects" in self.params[action_numbered]:
                        if isinstance(self.params[action_numbered]["redirects"], dict):
                            redirects = " \\\n\t".join(
                                [key + " " + (str(val) if val else "")
                                 for key, val
                                 in self.params[action_numbered]["redirects"].items()])
                        else:
                            redirects = self.params[action_numbered]["redirects"]
                else:
                    redirects = self.params[action_numbered]
                # Get action region
                if region:
                    if "region" in self.samtools_params[action]:
                        if self.samtools_params[action]["region"] is None:
                            self.write_warning(
                                "Region passed to tool '{action}' that does not accept regions. Removing!!".format(
                                    action=action))
                            region = ""

                        elif self.samtools_params[action]["region"] == "END":
                            pass
                        else:   # self.samtools_params[action]["region"] == "-r". for eaxmple
                            redirects = " \\\n\t".join([redirects,
                                                  "{arg} {region}".format(
                                                      arg=self.samtools_params[action]["region"],
                                                      region=region)])
                            region = ""

                    else:
                        self.write_warning("Region passed to tool '{action}' that does not accept regions. Removing!!".format(action=action))
                        region = ""

                if action in "flags split targetcut rmdup".split(" "):
                    raise AssertionExcept("Tool {action} is not supported by the module".format(action=action))

                if action in ["view","sort", "index", "flagstat","stats","idxstats","depth","fasta","fastq", "mpileup"]:
                    output_type = self.get_action_output_type(action, active_type, redirects)
                    # outfile_cmd = self.samtools_params_dict[action]["outfile"]
                    outfile_cmd = self.samtools_params[action]["outfile"]
                    if re.search("system", outfile_cmd):
                        raise Exception("WARNING: os.system call in outfile definition! BEWARE")
                    outfile = eval(outfile_cmd, globals(), locals())

                    cmd = self.samtools_params[action]["script"].format(action=action,
                                                                             env_path=self.get_script_env_path(),
                                                                             active_file=active_files[active_type],
                                                                             params="" if not redirects else "\n\t" + redirects + " \\",
                                                                             region="" if not region else "\\\n\t" + region,
                                                                             outfile=temp_use_dir + outfile)
                    self.script += """\
###########
# Running samtools {action}
#----------------
{cmd}
""".format(cmd=cmd,
           action=action)

                    _locals = dict()

                    # Get relevant local variables into _locals
                    for k in "action,active_type,output_type,active_files,files2keep,temp_use_dir,sample_dir," \
                             "outfile,sample".split(","):
                        _locals[k] = locals()[k]
                    active_type, active_files, files2keep = self.file_management(**_locals)

                else:
                    raise AssertionExcept("Action '{action}' not defined yet".format(action=action))

# TODO: incorporate filtering in this version of samtools
#             # The following can be merged into the main 'view' section
#             if "filter_by_tag" in list(self.params.keys()):
#                 # outfile = os.path.basename(active_file) + filter_suffix
#                 outfile = filter_suffix.join(os.path.splitext(os.path.basename(active_file)))
#
#                 self.script += """\
# ###########
# # Filtering BAM
# #----------------
#
# {env_path} view \\
# \t-h \\
# \t{active_file} | \\
# \tawk '$0 ~\"(^@)|({query})\"' | \\
# \t{env_path} view \\
# \t-bh \\
# \t-o {outfile} \\
# \t-
#
# {rm_unfilt}
# """.format(env_path=self.get_script_env_path(),
#            active_file=active_file,
#            query=self.params["filter_by_tag"],
#            outfile=temp_use_dir+outfile,
#            rm_unfilt="# Removing unfiltered BAM\nrm -rf "+active_file if "del_unfiltered" in list(self.params.keys()) else "")
#
#                 # Storing filtered and unfiltered bams:
#                 self.sample_data[sample]["unfiltered_bam"] = active_file
#                 self.sample_data[sample]["bam"] = sample_dir + outfile
#                 self.stamp_file(self.sample_data[sample]["bam"])
#                 active_file = temp_use_dir + outfile


            # files2keep = " \\\n\t".join(list(active_files.values()))
            files2keep = " \\\n\t".join(list(files2keep))
            self.script += """\
###########
# Copying final files to final location
#----------------
mv {files} \\
\t{dir}

rm -rf {tempdir}
""".format(files=files2keep,
           dir=use_dir,
           tempdir=temp_use_dir)


            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()



    def get_action_output_type(self, action, active_type, redirects):
        if action == "view":
            if re.search("\-\w*b", redirects):
                return "bam"
            elif re.search("\-\w*C", redirects):
                return "cram"
            else:
                return "sam"
        elif action == "sort":
            if re.search("\-\w*O", redirects):
                # TODO: get type from -O value
                pass
                # output_type = "bam"
            else:
                return "bam"
        elif action == "index":
            if active_type == "bam":
                return "bai"
            elif active_type == "cram":
                return "crai"
            else:
                raise AssertionExcept("No 'bam' or 'cram' for 'samtools index'", sample)
        elif action in ["flagstat", "stats", "idxstats", "depth"]:
            return action
        elif action in ["fasta", "fastq"]:
            return active_type
        elif action in ["mpileup"]:
            if re.search("\-\w*v", redirects) or re.search("\-\-VCF", redirects):
                return "vcf"
            elif re.search("\-\w*g", redirects) or re.search("\-\-BCF", redirects):
                return "bcf"
            else:
                return "mpileup"


    def file_management(self,action,active_type,output_type,active_files,files2keep,temp_use_dir,sample_dir,outfile,sample):

        if action in ["view", "sort", "index", "mpileup"]:
            active_files[output_type] = temp_use_dir + outfile
            self.sample_data[sample][output_type] = sample_dir + outfile
            if action in self.params["keep_output"]:
                files2keep.append(temp_use_dir + outfile)
            self.stamp_file(self.sample_data[sample][output_type])
            if action in ["view", "sort"]:
                active_type = output_type

        elif action in ["flagstat", "stats", "idxstats", "depth"]:
            active_files[output_type] = temp_use_dir + outfile
            self.sample_data[sample][active_type + "." + action] = sample_dir + outfile
            self.stamp_file(self.sample_data[sample][active_type + "." + action])
            if action in self.params["keep_output"]:
                files2keep.append(temp_use_dir + outfile)

        elif action in ["fasta","fastq"]:

            active_files[action + ".F"] = "%s%s.F.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
            active_files[action + ".R"] = "%s%s.R.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
            active_files[action + ".S"] = "%s%s.S.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
            # Storing and Stamping files
            self.sample_data[sample][action + ".F"] = "%s%s.F.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
            self.sample_data[sample][action + ".R"] = "%s%s.R.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
            self.sample_data[sample][action + ".S"] = "%s%s.S.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
            self.stamp_file(self.sample_data[sample][action + ".F"])
            self.stamp_file(self.sample_data[sample][action + ".R"])
            self.stamp_file(self.sample_data[sample][action + ".S"])
            if action in self.params["keep_output"]:
                files2keep.append(temp_use_dir + outfile)

        return active_type,active_files,files2keep




#
#
#     samtools_params_dict = \
#         {
#             "addreplacerg": None,
#             "bedcov": None,
#             "calmd": None,
#             "cat": None,
#             "collate": None,
#             "depad": None,
#             "depth": {
#                 "outfile": "\".\".join([os.path.basename(active_files[active_type]), output_type, \"txt\"])",
#                 "script": "{env_path}{action} \\{params}\n\t{active_file}  \\\n\t> {outfile}\n"
#             },
#             "dict": None,
#             "faidx": None,
#             "fasta": {
#                 "outfile": "os.path.basename(active_files[active_type])",
#                 "script": """\
# {env_path}{action} \\{params}
# \t-1  {outfile}.F.{action} \\
# \t-2  {outfile}.R.{action} \\
# \t-s  {outfile}.S.{action} \\
# \t{active_file}\n"""
#             },
#             "fastq": {
#                 "outfile": "os.path.basename(active_files[active_type])",
#                 "script": """\
# {env_path}{action} \\{params}
# \t-1  {outfile}.F.{action} \\
# \t-2  {outfile}.R.{action} \\
# \t-s  {outfile}.S.{action} \\
# \t{active_file}\n"""
#             },
#             "fixmate": None,
#             "flags": None,
#             "flagstat": {
#                 "outfile": "\".\".join([os.path.basename(active_files[active_type]), output_type, \"txt\"])",
#                 "script": "{env_path}{action} \\{params}\n\t{active_file}  \\\n\t> {outfile}\n"
#             },
#             "fqidx": None,
#             "idxstats": {
#                 "outfile": "\".\".join([os.path.basename(active_files[active_type]), output_type, \"txt\"])",
#                 "script": "{env_path}{action} \\{params}\n\t{active_file}  \\\n\t> {outfile}\n"
#             },
#             "index": {
#                 "outfile": "'{fn}.{ext}'.format(ext=output_type, fn=active_files[active_type])",
#                 "script": """\
# {env_path}{action} \\{params}
# \t{active_file}\n"""
#             },
#             "markdup": None,
#             "merge": None,
#             "mpileup": None,
#             "paired-end": None,
#             "phase": None,
#             "quickcheck": None,
#             "reheader": None,
#             "rmdup": None,
#             "sort": {
#                 "outfile": "os.path.splitext(os.path.basename(active_files[active_type]))[0] + \".sort.\" + output_type",
#                 "script": "{env_path}{action} \\{params}\n\t-o {outfile} \\\n\t{active_file}\n"
#             },
#             "split": None,
#             "stats": {
#                 "outfile": "\".\".join([os.path.basename(active_files[active_type]), output_type, \"txt\"])",
#                 "script": """\
# {env_path}{action} \\{params}
# \t{active_file}  \\
# \t> {outfile}"""
#             },
#             "targetcut": None,
#             "tview": None,
#             "view": {
#                 "input": {
#                     "flag": "stdin",
#                     "type": [
#                         "sam",
#                         "bam",
#                         "cram"
#                     ]
#                 },
#                 "outfile": "os.path.splitext(os.path.basename(active_files[active_type]))[0] + '.view.' + output_type",
#                 "script": """/
# {env_path}{action} \\{params}
# \t-o {outfile} \\
# \t{active_file} {region} """
#             }
#         }