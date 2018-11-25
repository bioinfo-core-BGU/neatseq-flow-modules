# -*- coding: UTF-8 -*-
"""
``samtools`` :sup:`*`
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

.. Note:: Order of samtools subprogram execution:

    The ``samtools`` programs are executed in the order above. It is up to you to have a sensible combination...


**Requires**:

* A SAM file in the following location:

    * ``sample_data[<sample>]["sam"]`` (for ``scope=sample``)
    * ``sample_data["sam"]`` (for ``scope=project``)

* Or a BAM file in:

    * ``sample_data[<sample>]["bam"]`` (for ``scope=sample``)
    * ``sample_data["bam"]`` (for ``scope=project``)

**Output**:

* Depending on the parameters, will put files in the following locations:

    * ``sample_data[<sample>]["bam"]``
    * ``sample_data[<sample>]["index"]``
    * ``sample_data[<sample>]["unfiltered_bam"]``
    * ``sample_data[<sample>]["unsorted_bam"]``
    * ``sample_data[<sample>]["flagstat"]``
    * ``sample_data[<sample>]["stats"]``
    * ``sample_data[<sample>]["idxstats"]``

* If ``fastq`` was called, will also create the following files:

    * ``self.sample_data[<sample>]["fastq.F"]``
    * ``self.sample_data[<sample>]["fastq.R"]``
    * ``self.sample_data[<sample>]["fastq.S"]``

* If ``fasta`` was called, will also create the following files:

    * ``self.sample_data[<sample>]["fasta.F"]``
    * ``self.sample_data[<sample>]["fasta.R"]``
    * ``self.sample_data[<sample>]["fasta.S"]``


.. Note:: If ``sample`` is set to ``project``, the above mentioned output files will be created in the project
    scope, e.g. ``sample_data["stats"]``..

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


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_samtools(Step):
       

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        if "type2use" in self.params:
            if not isinstance(self.params["type2use"],str) or self.params["type2use"] not in ["sam","bam"]:
                raise AssertionExcept("'type2use' must be either 'sam' or 'bam'")

        # Setting default scope to sample
        if "scope" not in self.params:
            self.params["scope"] = "sample"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Set list of samples to go over. Either self.sample_data["samples"] for sample scope
        # or ["project_data"] for project scope
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Check that a sam or bam exists
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
                raise AssertionExcept("Neither BAM nor SAM file exist.\n", sample)

        if "type2use" in self.params:
            if self.params["type2use"] not in self.sample_data[sample]:
                raise AssertionExcept("No file of type '{type}' exists.".format(type=self.params["type2use"]),
                                      sample)


    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

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

            active_file = self.sample_data[sample][self.file2use]

            filter_suffix = ".filt.bam"
            sort_suffix = ".srt.bam"
            index_suffix = ".bai"

            if "view" in self.params.keys():

                output_type = "bam" if re.search("\-\w*b", self.params["view"]) else "sam"
                outfile = ".".join([os.path.basename(active_file), output_type])

                self.script += """\
###########
# Running samtools view
#----------------

{env_path} view \\{params}
\t-o {outfile} \\
\t{active_file} 

""".format(env_path=self.get_script_env_path(),
                   active_file=active_file,
                   params="" if not self.params["view"] else "\n\t" + self.params["view"] + " \\",
                   outfile=use_dir + outfile)

                active_file = use_dir + outfile
                self.sample_data[sample][output_type] = sample_dir + outfile
                self.stamp_file(self.sample_data[sample][output_type])

            # If target of view is sam, terminating script. All others work on bam only.
                if output_type == "sam":
                    self.write_warning("""
                    Output from samtools view is SAM. Not proceeding further.
                    To produce a BAM, make sure to include the -b flag in the samtools view parameters.""")
                    # If sam output, can't proceed with rest of commands which require bam input_file:
                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir, sample_dir)
                    self.create_low_level_script()
                    continue
            else:
                # view not passed
                # If source is SAM, terminate with error
                if self.file2use == "sam":
                    raise AssertionExcept("Source file is 'sam', you must include 'view' in your oprations")
                # else, create local link to BAM and set active file accordingly
                self.script += """\
##########
# Making local link to original bam file:
#----------
cp -s {active_file} {here}

""".format(active_file=active_file,
           here=use_dir)
                active_file = sample_dir + os.path.basename(active_file)

            # The following can be merged into the main 'view' section
            if "filter_by_tag" in self.params.keys():
                outfile = os.path.basename(active_file) + filter_suffix

                self.script += """\
###########
# Filtering BAM
#----------------

{env_path} view \\
\t-h \\
\t{active_file} | \\  
\tawk '$0 ~\"(^@)|({query})\"' | \\
\t{env_path} view \\
\t-bh \\
\t-o {outfile} \\
\t- 

{rm_unfilt}
""".format(env_path=self.get_script_env_path(),
           active_file=active_file,
           query=self.params["filter_by_tag"],
           filt_suff=filter_suffix,
           outfile=use_dir+outfile,
           rm_unfilt="# Removing unfiltered BAM\nrm -rf "+active_file if "del_unfiltered" in self.params.keys() else "")

                # Storing filtered and unfiltered bams:
                self.sample_data[sample]["unfiltered_bam"] = active_file
                self.sample_data[sample]["bam"] = sample_dir + outfile
                self.stamp_file(self.sample_data[sample]["bam"])
                active_file = use_dir + outfile

            if "sort" in self.params.keys():
                if "bam" not in self.sample_data[sample]:
                    raise AssertionExcept("Can't run 'sort', as no BAM is defined", sample)
                outfile = os.path.basename(active_file) + sort_suffix

                self.script += """\
###########
# Sorting BAM
#----------------
{env_path}sort \\{params}
\t-o {outf} \\
\t{active_file}    
            
{rm_unsort}

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params["sort"] else "\n\t"+self.params["sort"]+" \\",
           outf=(use_dir + outfile),
           active_file=active_file,
           rm_unsort="# Removing unsorted BAM\nrm -rf "+active_file if "del_unsorted" in self.params.keys() else "")

                # Storing sorted bam in 'bam' slot and unsorted bam in unsorted_bam slot
                self.sample_data[sample]["unsorted_bam"] = active_file
                self.sample_data[sample]["bam"] = sample_dir + outfile
                self.stamp_file(self.sample_data[sample]["bam"])
                active_file = use_dir + outfile

            if "index" in self.params.keys():

                self.script += """\
###########
# Indexing BAM
#----------------
{env_path}index \\{params}
\t{active_file}    

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params["index"] else "\n\t" + self.params["index"] + " \\",
           active_file=active_file)

                self.sample_data[sample]["index"] = sample_dir + os.path.basename(active_file) + index_suffix
                self.stamp_file(self.sample_data[sample]["index"])


            for comm in ["flagstat","stats","idxstats"]:
                if comm in self.params.keys():
                    outfile = ".".join([os.path.basename(active_file), comm])

                    self.script += """\
###########
# Calculating {comm}
#----------------
{env_path}{comm} \\{params}
\t{active_file}  \\  
\t> {outfile}

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params[comm] else "\n\t" + self.params[comm] + " \\",
           active_file=active_file,
           comm=comm,
           outfile=use_dir+outfile)

                    self.sample_data[sample][comm] = sample_dir + outfile
                    self.stamp_file(self.sample_data[sample][comm])

            # Adding code for fastq or fasta extraction from bam:
            for type in (set(self.params.keys()) & set(["fasta","fastq"])):

                if "fastq.F" in self.sample_data[sample]:
                    readspart = """\
-1  {readsF} \\
\t-2  {readsR} \
""".format(readsF=(active_file + ".F." + type),
           readsR=(active_file + ".R." + type))
                else:
                    readspart = """\
-s  {readsS} \
""".format(readsS=(active_file + ".S." + type))

                # -0 and mixed paired-single not supported yet

                self.script += """\
###########
# Extracting fastq files from BAM:
#----------------
{env_path}{type} \\{params}
\t{readspart} \\
\t{active_file}

""".format(env_path=self.get_script_env_path(),
           params="" if not self.params[type] else "\n\t" + self.params[type] + " \\",
           readspart=readspart,
           type=type,
           active_file=active_file)

                # Storing and Stamping files
                if "fastq.F" in self.sample_data[sample]:
                    self.sample_data[sample][type+".F"] = "%s%s.F.%s" % (sample_dir, os.path.basename(active_file), type)
                    self.sample_data[sample][type+".R"] = "%s%s.R.%s" % (sample_dir, os.path.basename(active_file), type)
                    self.stamp_file(self.sample_data[sample][type+".F"])
                    self.stamp_file(self.sample_data[sample][type+".R"])
                else:
                    self.sample_data[sample][type+".S"] = "%s%s.S.%s" % (sample_dir, os.path.basename(active_file), type)
                    self.stamp_file(self.sample_data[sample][type+".S"])

            if "del_sam" in self.params.keys() and "sam" in self.sample_data[sample]:
                self.script += """\
###########
# Removing SAM
#----------------

rm -rf {sam}

""".format(sam=self.sample_data[sample]["sam"])

            self.local_finish(use_dir,sample_dir)
            self.create_low_level_script()

