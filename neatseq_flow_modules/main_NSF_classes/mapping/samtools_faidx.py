# -*- coding: UTF-8 -*-
"""
``samtools_faidx`` :sup:`*`
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

__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_samtools_faidx(Step):
       

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

            if "faidx" in self.params:

                fastatypes = list({"fasta.nucl", "fasta.prot"} & set(self.sample_data[sample]))
                if len(fastatypes) > 1:
                    if "fasta2use" not in self.params or self.params["fasta2use"] not in ["nucl", "prot"]:
                        raise AssertionExcept("You have both fasta.nucl and fasta.prot defined. Please select "
                                              "which to use by setting 'fasta2use' to 'nucl' or 'prot' ")
                    else:
                        try:
                            self.fasta2use = "fasta." + self.params["fasta2use"]
                        except KeyError:
                            raise AssertionExcept("You set fasta2use={f2use} but it does not "
                                                  "exist!".format(f2use=self.fasta2use), sample)
                elif len(fastatypes) == 0:
                    raise AssertionExcept("No fasta files exist!", sample)
                else:
                    print(fastatypes)
                    self.fasta2use = fastatypes[0]

            try:
                self.file2use
            except AttributeError:
                # If no sam/bam requiring actions exist, use fasta a file2use
                self.file2use = self.fasta2use




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
                if action not in self.samtools_params:
                    continue
                # print("------------------")
                # print(action)
                # print(active_files)

                if action in "flags split targetcut rmdup".split(" "):
                    raise AssertionExcept("Tool {action} is not supported by the module".format(action=action))


                if action == "view":

                    if re.search("\-\w*b", self.params[action]):
                        output_type = "bam"
                    elif re.search("\-\w*C", self.params[action]):
                        output_type = "cram"
                    else:
                        output_type = "sam"
                    # outfile = ".".join([os.path.basename(active_files[active_type]), output_type])
                    outfile = ("."+action+".").join([os.path.splitext(os.path.basename(active_files[active_type]))[0],output_type])

                    self.script += """\
###########
# Running samtools {action}
#----------------

{env_path} {action} \\{params}
\t-o {outfile} \\
\t{active_file} {region} 

""".format(action=action,
           env_path=self.get_script_env_path(),
           active_file=active_files[active_type],
           params="" if not self.params[action] else "\n\t" + self.params[action] + " \\",
           region="" if not "region" in self.params else "\\\n\t" + self.params["region"],
           outfile=temp_use_dir + outfile)

                    active_type = output_type
                    active_files[active_type] = temp_use_dir + outfile
                    self.sample_data[sample][output_type] = sample_dir + outfile
                    self.stamp_file(self.sample_data[sample][output_type])


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

            # if "sort" in list(self.params.keys()):

                if action == "sort":

                    if re.search("\-\w*O", self.params[action]):
                        # TODO: get type from -O value
                        pass
                        # output_type = "bam"
                    else:
                        output_type = "bam"
                    # outfile = ".".join([os.path.basename(active_files[active_type]), output_type])
                    # outfile = ("."+action).join(os.path.splitext(os.path.basename(active_files[active_type])))
                    outfile = ("."+action+".").join([os.path.splitext(os.path.basename(active_files[active_type]))[0],output_type])

                    # if "bam" not in self.sample_data[sample]:
                #     raise AssertionExcept("Can't run 'sort', as no BAM is defined", sample)
                # outfile = os.path.basename(active_file) + sort_suffix


                    self.script += """\
###########
# Running samtools {action}
#----------------
{env_path} {action} \\{params}
\t-o {outf} \\
\t{active_file}    
            
{rm_unsort}

""".format(action=action,
           env_path=self.get_script_env_path(),
           params="" if not self.params[action] else "\n\t"+self.params[action]+" \\",
           outf=(temp_use_dir + outfile),
           active_file=active_files[active_type],
           rm_unsort="# Removing unsorted BAM\nrm -rf " + active_files[active_type] if "del_unsorted" in list(self.params.keys()) else "")

                    # # Storing sorted bam in 'bam' slot and unsorted bam in unsorted_bam slot
                    # self.sample_data[sample]["unsorted_bam"] = active_file
                    # self.sample_data[sample]["bam"] = sample_dir + outfile
                    # self.stamp_file(self.sample_data[sample]["bam"])
                    # active_file = temp_use_dir + outfile
                    #
                    #

                    active_type = output_type
                    active_files[active_type] = temp_use_dir + outfile
                    self.sample_data[sample][output_type] = sample_dir + outfile
                    self.stamp_file(self.sample_data[sample][output_type])


                if action == "index":

                    if active_type=="bam":
                        output_type = "bai"
                    elif active_type == "cram":
                        output_type = "crai"
                    else:
                        raise AssertionExcept("No 'bam' or 'cram' for 'samtools index'", sample)
                    outfile = "{fn}.{ext}".format(ext=output_type, fn=active_files[active_type])


                    self.script += """\
###########
# Indexing BAM
#----------------
{env_path}{action} \\{params}
\t{active_file}    

""".format(action=action,
           env_path=self.get_script_env_path(),
           params="" if not self.params[action] else "\n\t" + self.params[action] + " \\",
           active_file=active_files[active_type])

                    # active_files[active_type] = temp_use_dir + outfile  - index does not change active type!
                    active_files[output_type] = outfile
                    self.sample_data[sample][output_type] = sample_dir + os.path.basename(outfile)
                    self.stamp_file(self.sample_data[sample][output_type])

                    # self.sample_data[sample]["bai"] = "{dir}{fn}.bai".format(dir=sample_dir,
                    #                                                          fn=os.path.basename(active_file))
                    # self.stamp_file(self.sample_data[sample]["bai"])


            # self.local_finish(temp_use_dir,sample_dir)
            # self.create_low_level_script()
            #
            # continue
                if action in ["flagstat","stats","idxstats","depth"]:
                    output_type = action
                    outfile = ".".join([os.path.basename(active_files[active_type]), output_type, "txt"])


                    self.script += """\
###########
# Calculating {action}
#----------------
{env_path}{action} \\{params}
\t{active_file}  \\
\t> {outfile}

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params[action] else "\n\t" + self.params[action] + " \\",
           active_file=active_files[active_type],
           action=action,
           outfile=temp_use_dir+outfile)

                    active_files[output_type] = temp_use_dir + outfile
                    self.sample_data[sample][active_type+"."+action] = sample_dir + outfile
                    self.stamp_file(self.sample_data[sample][active_type+"."+action])




                if action in ["fasta","fastq"]:

                    if active_type == "bam":
                        output_type = "bai"
                    elif active_type == "cram":
                        output_type = "crai"
                    else:
                        raise AssertionExcept("No 'bam' or 'cram' for 'samtools index'", sample)
                    outfile = ".".join([os.path.basename(active_files[active_type]), output_type])

            # # Adding code for fastq or fasta extraction from bam:
            # for type in (set(self.params.keys()) & set(["fasta","fastq"])):

#                     if "fastq.F" in self.sample_data[sample]:
#                         readspart = """\
# -1  {readsF} \\
# \t-2  {readsR} \
# """.format(readsF=(active_files[active_type] + ".F." + type),
#            readsR=(active_files[active_type] + ".R." + type))
#                     else:
#                         readspart = """\
# -s  {readsS} \
# """.format(readsS=(active_files[active_type] + ".S." + type))
                    readspart = """\
-1  {readsF} \\
\t-2  {readsR} \\
\t-s  {readsS} \
""".format(readsS=(active_files[active_type] + ".S." + action),
           readsF=(active_files[active_type] + ".F." + action),
           readsR=(active_files[active_type] + ".R." + action))


                    self.script += """\
###########
# Extracting fastq files from BAM:
#----------------
{env_path}{action} \\{params}
\t{readspart} \\
\t{active_file}

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params[action] else "\n\t" + self.params[action] + " \\",
           readspart=readspart,
           action=action,
           active_file=active_files[active_type])

                    active_files[action+".F"] = "%s%s.F.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
                    active_files[action+".R"] = "%s%s.R.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
                    active_files[action+".S"] = "%s%s.S.%s" % (temp_use_dir, os.path.basename(active_files[active_type]), action)
                    # Storing and Stamping files
                    self.sample_data[sample][action+".F"] = "%s%s.F.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
                    self.sample_data[sample][action+".R"] = "%s%s.R.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
                    self.sample_data[sample][action+".S"] = "%s%s.S.%s" % (sample_dir, os.path.basename(active_files[active_type]), action)
                    self.stamp_file(self.sample_data[sample][action+".F"])
                    self.stamp_file(self.sample_data[sample][action+".R"])
                    self.stamp_file(self.sample_data[sample][action+".S"])

                # --------------------------------------- faidx
                if action == "faidx":

                    fastatypes = list({"fasta.nucl", "fasta.prot"} & set(self.sample_data[sample]))
                    if len(fastatypes)>1:
                        if "fasta2use" not in self.params or self.params["fasta2use"] not in ["nucl","prot"]:
                            raise AssertionExcept("You have both fasta.nucl and fasta.prot defined. Please select "
                                                  "which to use by setting 'fasta2use' to 'nucl' or 'prot' ")
                        else:
                            try:
                                fasta2use = "fasta." + self.params["fasta2use"]
                                active_files[fasta2use] = self.sample_data[sample][fasta2use]
                            except KeyError:
                                raise AssertionExcept("You set fasta2use={f2use} but it does not "
                                                      "exist!".format(f2use=fasta2use),sample)
                    elif len(fastatypes) ==0:
                        raise AssertionExcept("No fasta files exist!", sample)
                    else:
                        print(fastatypes)
                        fasta2use = fastatypes[0]
                        active_files[fasta2use] = self.sample_data[sample][fasta2use]

                    output_type = "fai"
                    outfile = "{fn}.{ext}".format(ext=output_type, fn=active_files[active_type])

                    self.script += """\
###########
# Indexing fasta
#----------------
cp -sf {src} {trg}
{env_path}{action} \\{params}
\t{active_file}    

""".format(src=active_files[active_type],
           trg=temp_use_dir,
           action=action,
           env_path=self.get_script_env_path(),
           params="" if not self.params[action] else "\n\t" + self.params[action] + " \\",
           active_file= temp_use_dir + os.path.basename(active_files[active_type]))

                    # active_files[active_type] = temp_use_dir + outfile  - index does not change active type!
                    active_files[output_type] = outfile
                    self.sample_data[sample][output_type] = sample_dir + os.path.basename(outfile)
                    self.stamp_file(self.sample_data[sample][output_type])

                    # self.sample_data[sample]["bai"] = "{dir}{fn}.bai".format(dir=sample_dir,
                    #                                                          fn=os.path.basename(active_file))
                    # self.stamp_file(self.sample_data[sample]["bai"])

            self.script += """\
###########
# Copying final files to final location
#----------------
mv {files} \\
\t{dir}

rm -rf {tempdir}
""".format(files=" \\\n\t".join(list(active_files.values())),
           dir=use_dir,
           tempdir=temp_use_dir)


            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()

            continue
#
#             if "del_sam" in list(self.params.keys()) and "sam" in self.sample_data[sample]:
#                 self.script += """\
# ###########
# # Removing SAM
# #----------------
#
# rm -rf {sam}
#
# """.format(sam=self.sample_data[sample]["sam"])
#
#             self.local_finish(use_dir,sample_dir)
#             self.create_low_level_script()




