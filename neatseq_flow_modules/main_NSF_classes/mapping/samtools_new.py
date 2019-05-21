# -*- coding: UTF-8 -*-
"""
``samtools_new`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for executing samtools on a SAM or BAM file.

.. Warning:: This module is in beta stage. Please report issues and we'll try solving them

.. attention:: The module was tested on samtools 1.9

Currently, the samtools programs included in the module are the following:

* ``view``
* ``sort``
* ``index``
* ``flagstat``
* ``stats``
* ``idxstats``
* ``depth``
* ``fastq/a``
* ``merge``
* ``mpileup``

.. Note:: Order of samtools subprogram execution:

    The ``samtools`` programs are executed in the order given in the parameter file!

    File types are passed from one program to the next based on the output type of the program.

    In order to execute one program more than once, append digits to the program name, *e.g.* ``sort2``, ``index3`` etc.

Arguments can be passed to the tools following the program name in the parameter file, *e.g.*::

    sort: -n -@ 10

Please do NOT pass input and output arguments - they are set by the module.

If you want to limit the program to a specific region, pass the program name a block with a 'region' section. If you want to set the region and pass some redirects, add a 'redirect' section as well. For example::

    mpileup:
        redirects:      --max-depth INT -v
        region:         chr2:212121-32323232

Some of the tools are defined only when the ``scope`` is ``sample``, namely ``merge`` and ``mpileup``. ``merge`` merges
the sample-wise BAM files into a project BAM file. ``mpileup`` creates a project VCF/BCF/mpileup file from the sample BAM files.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* A SAM file in the following location:

    * ``sample_data[<sample>]["sam"]`` (for ``scope=sample``)
    * ``sample_data["project_data"]["sam"]`` (for ``scope=project``)

* Or a BAM file in:

    * ``sample_data[<sample>]["bam"]`` (for ``scope=sample``)
    * ``sample_data["project_data"]["bam"]`` (for ``scope=project``)

.. Note:: If both ``BAM`` and ``SAM`` files exist, select the one to use with ``type2use`` (see section *Parameters that can be set*).

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Depending on the parameters, will put files in different types (*e.g.* ``bam``, ``cram``, ``sam``, ``bam``, ``bai``, ``crai``, ``vcf``, ``bcf``, ``mpileup``, ``fasta.{F,R,S}``, ``fastq.{F,R,S}``)
Please use ``stop_and_show`` to see the types produced by your instance of ``samtools_new``.

.. Note:: If ``scope`` is set to ``project``, the above mentioned output files will be created in the project scope.

.. Note:: ``merge`` and ``mpileup`` are only defined when ``scope`` is ``sample``. See above

By default, all files are saved. To keep only the output from specific programs, add a ``keep_output`` section containing a list of programs for which the output should be saved. **All other files will be discarded**.
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
    "region", "", "A region to limit the region-limitable programs, such as ``view``, ``merge``, ``mpileup``, etc.."
    "type2use","sam|bam","Type of file to use. Must exist in scope"
    "keep_output", "[sort, view, sort2]", "A list of programs for which to store the output files. By deafult, all files are saved.
..    "filter_by_tag", "*e.g.*: NM:i:[01]", "Filter BAM by one of the tags. Use an awk-compliant regular expression. In this example, keep only lines where the edit distance is 0 or 1. This is an experimental feature and should be used with caution..."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sam_bwt1:
        module:             samtools_new
        base:               bwt1
        script_path:        {Vars.paths.samtools}
        qsub_params:
            -pe:            shared 20
        region:             chr2:212121-32323232
        scope:              sample
        # First 'view'. Use FLAG to filter alignments:
        view:               -uh  -q 30 -@ 20 -F 4 -O bam
        # First 'sort'. Sort by coordinates:
        sort:               -@ 20
        # Second 'view'. Use region to filter alignments:
        view2:
            redirects:      -buh  -q 30 -@ 20
            region:         chr2:212121-32323232
        index:
        flagstat:
        stats:              --remove-dups
        idxstats:
        # Second 'sort'. Sort by name:
        sort2:               -n -@ 20
        # Get sequences from name-sorted BAM file:
        fastq:
        # Merge BAM name sorted BAM files
        merge:
            region:         chr2:212121-32323232
        # Create VCF from Merge BAM name sorted BAM files
        mpileup:
            redirects:      --max-depth INT -v
            region:         chr2:212121-32323232
        keep_output:        [sort, view, index, flagstat, stats, fastq, mpileup, merge]
        # stop_and_show:


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

        not_supported_actions = "dict faidx fqidx flags split targetcut rmdup".split(" ")
        not_supported_actions = list(set(self.params) & set(not_supported_actions))
        if len(not_supported_actions) > 0:
            raise AssertionExcept("These tools are not supported by this module: {tools}".format(tools=", ".join(not_supported_actions)))

        self.supported_actions = "view sort index flagstat stats idxstats depth fasta fastq".split(" ")
        self.project_actions = "merge mpileup".split(" ")


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

        if any([action in self.project_actions for action in self.params]) and self.params["scope"] != "sample":
            raise AssertionExcept("{tools} are defined only for sample-scope instances.".
                                  format(tools=", ".join([action in self.project_actions for action in self.params])))

        if "keep_output" not in self.params:
            self.params["keep_output"] = list(set(self.params.keys()) & set(self.samtools_params.keys()))
        # if "del_output" not in self.params:
        #     self.params["del_output"] = list()

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
        # print(self.active_type)
        # sys.exit()
        # -------------- samtools merge ------------------
        sample = "project_data"
        self.script = ""

        files2keep = list()

        for action in self.params:
            # This is to enable passing each tool more than once:
            # Repeated actions have digits appended to them
            action_numbered = action
            action = re.sub("\d+$", "", action)

            if action not in self.samtools_params:
                continue

            if action in self.project_actions:
                # Get action redirects (if exist)
                redirects, region = self.region_and_redirects(action, action_numbered)
                # print("{name} --- {redirs}".format(name=action, redirs=redirects))

                if not self.active_type == "bam":
                    raise AssertionExcept("samtools {action} requires a bam file, but the active file type is "
                                          "{active_type}. Please check the arguments to previous samtools commands!".
                                          format(active_type=self.active_type,
                                                 action=action))
                redirects, region = self.region_and_redirects(action, action_numbered)

                active_files = dict(zip((self.active_type,),
                                        (" ".join([self.sample_data[sample][self.active_type] for sample in
                                                  self.sample_data["samples"]]),)))

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                sample_dir = self.make_folder_for_sample(sample)
                use_dir = self.local_start(self.base_dir)
                # outfile = self.sample_data["Title"] + ".merged.bam"
                output_type = self.get_action_output_type(action, redirects)
                outfile_cmd = self.samtools_params[action]["outfile"]
                outfile = eval(outfile_cmd, globals(), locals())


                cmd = self.samtools_params[action]["script"].format(action=action,
                                                                    env_path=self.get_script_env_path(),
                                                                    active_file=active_files[self.active_type],
                                                                    params="" if not redirects else "\n\t" + redirects + " \\",
                                                                    region="" if not region else "\\\n\t" + region,
                                                                    outfile=use_dir + outfile)
                self.script += """\
###########
# Running samtools {action}
#----------------
{cmd}
""".format(cmd=cmd,
           action=action)

                _locals = dict()

                # Get relevant local variables into _locals
                for k in "action,action_numbered,output_type,active_files,files2keep,use_dir," \
                         "sample_dir,outfile,sample".split(","):
                    _locals[k] = locals()[k]
                self.active_type, active_files, files2keep = self.file_management(**_locals)


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

            # active_file = self.sample_data[sample][self.file2use]
            active_files = dict(zip((self.file2use,),(self.sample_data[sample][self.file2use],)))
            self.active_type = self.file2use

            files2keep = list()

            # Starting off with local link to active file:
            self.script += """\
##########
# Making local link to original bam file: (-f to force)
#----------
cp -fs \\
\t{active_file} \\
\t{use_dir}

""".format(active_file=active_files[self.active_type],
           use_dir=use_dir)

            active_files[self.active_type] = use_dir + os.path.basename(active_files[self.active_type])
            self.sample_data[sample][self.active_type] = sample_dir + os.path.basename(active_files[self.active_type])

            for action in self.params:
                # This is to enable passing each tool more than once:
                # Repeated actions have digits appended to them
                action_numbered = action
                action = re.sub("\d+$","",action)

                if action not in self.samtools_params:
                    continue

                # Get action redirects (if exist)
                redirects, region = self.region_and_redirects(action, action_numbered)


                if action in self.supported_actions:
                    output_type = self.get_action_output_type(action, redirects)
                    # outfile_cmd = self.samtools_params_dict[action]["outfile"]
                    outfile_cmd = self.samtools_params[action]["outfile"]
                    if re.search("system", outfile_cmd):
                        raise Exception("WARNING: os.system call in outfile definition! BEWARE")
                    outfile = eval(outfile_cmd, globals(), locals())

                    cmd = self.samtools_params[action]["script"].format(action=action,
                                                                             env_path=self.get_script_env_path(),
                                                                             active_file=active_files[self.active_type],
                                                                             params="" if not redirects else "\n\t" + redirects + " \\",
                                                                             region="" if not region else "\\\n\t" + region,
                                                                             outfile=use_dir + outfile)
                    self.script += """\
###########
# Running samtools {action}
#----------------
{cmd}
""".format(cmd=cmd,
           action=action)

                    _locals = dict()

                    # Get relevant local variables into _locals
                    for k in "action,action_numbered,output_type,active_files,files2keep,use_dir," \
                             "sample_dir,outfile,sample".split(","):
                        _locals[k] = locals()[k]
                    self.active_type, active_files, files2keep = self.file_management(**_locals)

                elif action in self.project_actions:
                    pass
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
            if files2keep:
                self.script += """\
###########
# Copying final files to final location
#----------------
mkdir {temp}
mv {files} \\
\t{temp}

find {usedir} -maxdepth 1 -mindepth 1 -not -name "temp" -exec rm -rf {{}} \;

mv {temp}/* \\
\t{usedir}

rmdir {temp}
""".format(files=files2keep,
           dir=use_dir,
           usedir=use_dir,
           temp=temp_use_dir)


            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()

    def region_and_redirects(self, action, action_numbered):
        redirects = ""
        region = self.params["region"] if "region" in self.params else ""  # Enable globally defined region
        if isinstance(self.params[action_numbered], dict):
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
                    pass  # Region is dealt with later
                else:  # self.samtools_params[action]["region"] == "-r". for eaxmple
                    if redirects:
                        redirects = " \\\n\t".join([redirects,
                                                    "{arg} {region}".format(
                                                        arg=self.samtools_params[action]["region"],
                                                        region=region)])
                    else:
                        redirects = "{arg} {region}".format(
                                                        arg=self.samtools_params[action]["region"],
                                                        region=region)
                    region = ""  # Region is added to redirects and emptied

            else:
                if "region" not in self.params:  # Issue the comment only if region is set per this action
                    self.write_warning(
                        "Region passed to tool '{action}' that does not accept regions. Removing!!".format(
                            action=action))
                region = ""
        return redirects, region

    def get_action_output_type(self, action, redirects):

        if redirects is not None:
            if re.search("\-\w*O\s", redirects):
                type = re.search("\-\w*O\s+(\w+)", redirects)
                if type.group(1).lower() not in ["bam", "sam", "cram"]:
                    raise AssertionExcept("Bad value for output format ({type})".format(type=type.group(1)))
                else:
                    return type.group(1).lower()

        if action == "view":
            if re.search("\-\w*b", redirects):
                return "bam"
            elif re.search("\-\w*C", redirects):
                return "cram"
            else:
                return "sam"
        elif action == "sort":
            return "bam"
        elif action == "index":
            if self.active_type == "bam":
                return "bai"
            elif self.active_type == "cram":
                return "crai"
            else:
                raise AssertionExcept("No 'bam' or 'cram' for 'samtools index'", sample)
        elif action in ["flagstat", "stats", "idxstats", "depth"]:
            return action
        elif action in ["fasta", "fastq"]:
            return self.active_type
        elif action in ["mpileup"]:
            if re.search("\-\w*v", redirects) or re.search("\-\-VCF", redirects):
                return "vcf"
            elif re.search("\-\w*g", redirects) or re.search("\-\-BCF", redirects):
                return "bcf"
            else:
                return "mpileup"
        elif action in ["merge"]:
            return "bam"

    def file_management(self,action,action_numbered,output_type,active_files,files2keep,use_dir,sample_dir,outfile,sample):

        if action in ["view", "sort", "index", "mpileup", "merge"]:
            active_files[output_type] = use_dir + outfile
            self.sample_data[sample][output_type] = sample_dir + outfile
            if action_numbered in self.params["keep_output"]:
                files2keep.append(use_dir + outfile)
                self.stamp_file(self.sample_data[sample][output_type])
            if action in ["view", "sort"]:
                self.active_type = output_type

        elif action in ["flagstat", "stats", "idxstats", "depth"]:
            active_files[output_type] = use_dir + outfile
            self.sample_data[sample][self.active_type + "." + action] = sample_dir + outfile
            if action_numbered in self.params["keep_output"]:
                self.stamp_file(self.sample_data[sample][self.active_type + "." + action])
                files2keep.append(use_dir + outfile)

        elif action in ["fasta","fastq"]:

            active_files[action + ".F"] = "%s%s.F.%s" % (use_dir, os.path.basename(active_files[self.active_type]), action)
            active_files[action + ".R"] = "%s%s.R.%s" % (use_dir, os.path.basename(active_files[self.active_type]), action)
            active_files[action + ".S"] = "%s%s.S.%s" % (use_dir, os.path.basename(active_files[self.active_type]), action)
            # Storing and Stamping files
            self.sample_data[sample][action + ".F"] = "%s%s.F.%s" % (sample_dir, os.path.basename(active_files[self.active_type]), action)
            self.sample_data[sample][action + ".R"] = "%s%s.R.%s" % (sample_dir, os.path.basename(active_files[self.active_type]), action)
            self.sample_data[sample][action + ".S"] = "%s%s.S.%s" % (sample_dir, os.path.basename(active_files[self.active_type]), action)
            self.stamp_file(self.sample_data[sample][action + ".F"])
            self.stamp_file(self.sample_data[sample][action + ".R"])
            self.stamp_file(self.sample_data[sample][action + ".S"])
            if action_numbered in self.params["keep_output"]:
                files2keep.append(active_files[action + ".F"])
                files2keep.append(active_files[action + ".R"])
                files2keep.append(active_files[action + ".S"])

        else:
            print("Unrecognised action '{action}'. Check!!!".formar(action=action))

        return self.active_type,active_files,files2keep


