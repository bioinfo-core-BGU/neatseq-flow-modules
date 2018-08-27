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
                
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if "scope" not in self.params:
            self.write_warning("No scope specified. Setting to 'sample'")
            self.params["scope"] = "sample"

        if self.params["scope"] == "project":
            self.step_sample_initiation_byproject()
        elif self.params["scope"] == "sample":
            self.step_sample_initiation_bysample()
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

    def step_sample_initiation_bysample(self):

        for sample in self.sample_data["samples"]:  # Getting list of samples out of samples_hash

            # Check that a sam or bam exists
            if "bam" in self.sample_data[sample] and "sam" in self.sample_data[sample]:
                if "type2use" in self.params:
                    self.file2use = self.params["type2use"]
                else:
                    raise AssertionExcept(
                        "Both BAM and SAM file types exist for sample. Specify which one to use with 'type2use'.\n",
                        sample)
            elif "bam" in self.sample_data[sample]:
                self.file2use = "bam"
            elif "sam" in self.sample_data[sample]:
                self.file2use = "sam"
            else:
                raise AssertionExcept("Neither BAM nor SAM file exist for sample.\n", sample)

        if "type2use" in self.params:
            if self.params["type2use"] not in self.sample_data[sample]:
                raise AssertionExcept("No file of type '{type}' exists for sample".format(type=self.params["type2use"]),
                                      sample)

        pass

    def step_sample_initiation_byproject(self):

        # Check that a sam or bam exists
        if "bam" in self.sample_data and "sam" in self.sample_data:
            if "type2use" in self.params:
                self.file2use = self.params["type2use"]
            else:
                raise AssertionExcept("Both BAM and SAM file types exist for project. "
                                      "Specify which one to use with 'type2use'.\n")
        elif "bam" in self.sample_data:
            self.file2use = "bam"
        elif "sam" in self.sample_data:
            self.file2use = "sam"
        else:
            raise AssertionExcept("Neither BAM nor SAM file exist for project.\n")

        if "type2use" in self.params:
            if self.params["type2use"] not in self.sample_data:
                raise AssertionExcept("No file of type '{type}' exists for project".
                                      format(type=self.params["type2use"]))

        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function

        """

        if self.params["scope"] == "project":
            self.build_scripts_byproject()
        elif self.params["scope"] == "sample":
            self.build_scripts_bysample()
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

    def build_scripts_bysample(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        # Each iteration must define the following class variables:
            # self.spec_script_name
            # self.script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            

            # sam_file = self.sample_data[sample]["sam"]
            input_file = self.sample_data[sample][self.file2use]
            # bam_name = self.sample_data[sample]["sam"] + ".bam"
            bam_name = os.path.basename(input_file) + ".bam"
            output_sam_name = os.path.basename(input_file) + ".sam" #might be used...
            
            active_file = input_file
            
            if "filter_by_tag" in self.params.keys():
                filtered_name = bam_name + ".filt.bam"
                sort_name = filtered_name + ".srt.bam"
            else:
                sort_name = bam_name + ".srt.bam"
            index_name = sort_name + ".bai"
            
            filter_suffix = ".filt.bam"
            sort_suffix = ".srt.bam"
            index_suffix = ".bai"

            
            if "view" in self.params.keys():
                self.script += "###########\n# Running samtools view:\n#----------------\n"
                self.script += "%s view \\\n\t" % self.get_script_env_path()
                if self.params["view"]:
                    self.script += "%s \\\n\t" % self.params["view"]

                tobam = re.search("\-\w*b",self.params["view"])
                if tobam:
                    self.script += "-o {outfile} \\\n\t{infile}\n\n".format(outfile=use_dir + bam_name,
                                                                             infile=active_file)
                    active_file = use_dir + bam_name
                    self.sample_data[sample]["bam"] = sample_dir + bam_name
                    self.stamp_file(self.sample_data[sample]["bam"])
                    
                else:
                    self.script += "-o {outfile} \\\n\t{infile}\n\n".format(outfile=use_dir + output_sam_name,
                                                                             infile=active_file)
                    self.sample_data[sample]["sam"] = sample_dir + output_sam_name
                    self.stamp_file(self.sample_data[sample]["sam"])

                    self.write_warning("""
Output from samtools view is SAM. Not proceeding further.
To produce a BAM, make sure to include the -b flag in the samtools view parameters.""")
                    # If sam output, can't proceed with rest of commands which require bam input_file:
                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                    self.create_low_level_script()
                    continue

            # The following can be merged into the main 'view' section
            if "filter_by_tag" in self.params.keys():

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
\t-o {active_file}{filt_suff} \\
\t- 

""".format(env_path=self.get_script_env_path(),
           active_file=active_file,
           query=self.params["filter_by_tag"],
           filt_suff=filter_suffix)

                # self.script += "###########\n# Filtering BAM\n#----------------\n"
                # self.script += "\n\n"
                # self.script += "%s view \\\n\t" % self.get_script_env_path()
                # self.script += "-h \\\n\t"
                # self.script += "%s | \\\n\t" % active_file    # self.sample_data[sample]["bam"]
                # self.script += "awk '$0 ~\"(^@)|(%s)\"' | \\\n\t" % self.params["filter_by_tag"]
                # self.script += "%s view \\\n\t" % self.get_script_env_path()
                # self.script += "-bh \\\n\t"
                # self.script += "-o %s \\\n\t" % (active_file + filter_suffix)   # use_dir + filtered_name)
                # self.script += "- \n\n"

                # If user requires than unsorted bam be removed:
                if "del_unfiltered" in self.params.keys():
                    self.script += "###########\n# Removing unfiltered BAM\n#----------------\n"
                    self.script += "\n\nrm -rf %s\n\n" % (active_file)

                # Storing filtered and unfiltered bams:
                self.sample_data[sample]["unfiltered_bam"] = sample_dir + os.path.basename(active_file)   #bam_name
                active_file = active_file + filter_suffix
                self.sample_data[sample]["bam"] = sample_dir + os.path.basename(active_file)
                self.stamp_file(self.sample_data[sample]["bam"])

                # The following is so that sort will work on the filtered file without playing around with the sort code:
                bam_name = filtered_name
                
            if "sort" in self.params.keys():
                # This permits running only sort and index, in case a bam file was produced in a differnet step.
                # if "view" in self.params.keys():
                    # bam_name = use_dir + bam_name
                # else:
                    # if "bam" in self.sample_data[sample].keys():
                        # bam_name = self.sample_data[sample]["bam"]
                    # elif "sam" in self.sample_data[sample].keys():
                        # bam_name = self.sample_data[sample]["sam"]
                        # self.write_warning("Can't find BAM but found SAM for sample. Using it instead of a BAM.\n", sample)
                    # else:
                        # raise AssertionExcept("Can't run sort without BAM file. Either include 'view' or use other BAM creating steps.\n",sample)
                self.script += """\
###########
# Sorting BAM
#----------------
{env_path}sort \\{params} 
\t-o {outf} \\
\t{inf}                
                """.format(env_path=self.get_script_env_path(),
                           params="" if not self.params["sort"] else "\n\t"+self.params["sort"]+" \\",
                           outf=(use_dir + os.path.basename(active_file) + sort_suffix),
                           inf=active_file)

                # If user requires than unsorted bam be removed:
                if "del_unsorted" in self.params.keys():
                    self.script += "###########\n# Removing unsorted BAM\n#----------------\n"
                    self.script += "\n\nrm -rf %s\n\n" % (active_file)

                # Storing sorted bam in 'bam' slot and unsorted bam in unsorted_bam slot
                self.sample_data[sample]["unsorted_bam"] = sample_dir + os.path.basename(active_file)
                active_file = active_file + sort_suffix
                self.sample_data[sample]["bam"] = sample_dir + os.path.basename(active_file)
                self.stamp_file(self.sample_data[sample]["bam"])

            if "index" in self.params.keys():
                self.script += "###########\n# Indexing BAM\n#----------------\n"
                self.script += "%s index \\\n\t" % self.get_script_env_path()
                if self.params["index"]:
                    self.script += "%s \\\n\t" % self.params["index"]
                self.script += "%s\n\n" % active_file   #(use_dir + bam_name)
                self.sample_data[sample]["index"] = sample_dir + os.path.basename(active_file) + index_suffix
                self.stamp_file(self.sample_data[sample]["index"])

            if "flagstat" in self.params.keys():
                self.script += "###########\n# Calculating BAM statistics:\n#----------------\n"
                self.script += "%s flagstat \\\n\t" % self.get_script_env_path()
                self.script += "%s \\\n\t" % active_file   #(use_dir + bam_name)
                self.script += "> %s.flagstat \n\n" % active_file    #(use_dir + bam_name)
                self.sample_data[sample]["flagstat"] = "%s%s.flagstat" % (sample_dir, os.path.basename(active_file))  #bam_name)
                self.stamp_file(self.sample_data[sample]["flagstat"])


            if "stats" in self.params.keys():
                self.script += "###########\n# Calculating BAM statistics:\n#----------------\n"
                self.script += "%s stats \\\n\t" % self.get_script_env_path()
                if self.params["stats"]:  # Adding parameters the user might pass
                    self.script += "%s \\\n\t" % self.params["stats"]
                self.script += "%s \\\n\t" % active_file  #(use_dir + bam_name)
                self.script += "> %s.stats \n\n" % active_file  #(use_dir + bam_name)
                self.sample_data[sample]["stats"] = "%s%s.stats" % (sample_dir, os.path.basename(active_file))
                self.stamp_file(self.sample_data[sample]["stats"])

            if "idxstats" in self.params.keys() and "idxstats" in self.sample_data[sample]:
                self.script += "###########\n# Calculating index statistics (idxstats):\n#----------------\n"
                self.script += "%s idxstats \\\n\t" % self.get_script_env_path()
                # idxstats has no uder defined parameters...
                self.script += "%s \\\n\t" % active_file  #(use_dir + bam_name)
                self.script += "> %s.idxstat.tab \n\n" % active_file  #(use_dir + bam_name)
                self.sample_data[sample]["idxstats"] = "%s%s.idxstat.tab" % (sample_dir, os.path.basename(active_file))  #bam_name)
                self.stamp_file(self.sample_data[sample]["idxstats"])

            # Adding code for fastq or fasta extraction from bam:
            for type in (set(self.params.keys()) & set(["fasta","fastq"])):
                self.script += """###########\n# Extracting fastq files from BAM:\n#----------------\n"""
                self.script += self.get_script_env_path()
                self.script += "{command} {params} \\\n\t".format(command = type, params = self.params[type])
                if "fastq.F" in self.sample_data[sample]:
                    self.script += "-1  {readsF} \\\n\t".format(readsF = (active_file + ".F." + type))
                    self.script += "-2  {readsR} \\\n\t".format(readsR = (active_file + ".R." + type))
                else:
                    self.script += "-0  {readsS} \\\n\t".format(readsS = (active_file + ".S." + type))
                # -s and mixed paired-single not supported yet
                self.script += "{bam} \n\n".format(bam = (active_file))
                
                # Storing and Stamping files
                if "fastq.F" in self.sample_data[sample]:
                    self.sample_data[sample]["%s.F"    % type] = "%s%s.F.%s" % (sample_dir, os.path.basename(active_file)   , type)   
                    self.sample_data[sample]["%s.R"    % type] = "%s%s.R.%s" % (sample_dir, os.path.basename(active_file)   , type)   
                    self.stamp_file(self.sample_data[sample]["%s.F"    % type] )
                    self.stamp_file(self.sample_data[sample]["%s.R"    % type] )
                else:
                    self.sample_data[sample]["%s.S"    % type] = "%s%s.S.%s" % (sample_dir, os.path.basename(active_file)   , type)   
                    self.stamp_file(self.sample_data[sample]["%s.S"    % type] )
                



                
            if "del_sam" in self.params.keys() and "sam" in self.sample_data[sample]:
                self.script += "###########\n# Removing SAM\n#----------------\n\n"
                self.script += "rm -rf %s\n\n" % self.sample_data[sample]["sam"]

            

            # self.stamp_dir_files(sample_dir)
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            
            self.create_low_level_script()

    def build_scripts_byproject(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        input_file = self.sample_data[self.file2use]
        bam_name = os.path.basename(input_file) + ".bam"
        output_sam_name = os.path.basename(input_file) + ".sam"  # might be used...

        active_file = input_file

        if "filter_by_tag" in self.params.keys():
            filtered_name = bam_name + ".filt.bam"
            sort_name = filtered_name + ".srt.bam"
        else:
            sort_name = bam_name + ".srt.bam"
        index_name = sort_name + ".bai"

        filter_suffix = ".filt.bam"
        sort_suffix = ".srt.bam"
        index_suffix = ".bai"

        if "view" in self.params.keys():
            self.script += "###########\n# Running samtools view:\n#----------------\n"
            self.script += "%s view \\\n\t" % self.get_script_env_path()
            if self.params["view"]:
                self.script += "%s \\\n\t" % self.params["view"]

            tobam = re.search("\-\w*b", self.params["view"])
            if tobam:
                self.script += "-o {outfile} \\\n\t{infile}\n\n".format(outfile=use_dir + bam_name,
                                                                        infile=active_file)
                active_file = use_dir + bam_name
                self.sample_data["bam"] = self.base_dir + bam_name
                self.stamp_file(self.sample_data["bam"])

            else:
                self.script += "-o {outfile} \\\n\t{infile}\n\n".format(outfile=use_dir + output_sam_name,
                                                                        infile=active_file)
                self.sample_data["sam"] = sample_dir + output_sam_name
                self.stamp_file(self.sample_data["sam"])

                self.write_warning(
                    "Output from samtools view is SAM. Not proceeding further.\n"
                    "To produce a BAM, make sure to include the -b flag in the samtools view parameters.\n")
                # If sam output, can't proceed with rest of commands which require bam input_file:
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)  # Sees to copying local files to final destination (and other stuff)
                self.create_low_level_script()
                return

        # The following can be merged into the main 'view' section
        if "filter_by_tag" in self.params.keys():

            self.script += """\
###########
# Filtering BAM
#----------------

{env_path} view \\
    -h \\
    {active_file} | \\  
    awk '$0 ~\"(^@)|({query})\"' | \\
    {env_path} view \\
    -bh \\
    -o {active_file}{filt_suff} \\
    - 
    
""".format(env_path=self.get_script_env_path(),
           active_file=active_file,
           query=self.params["filter_by_tag"],
           filt_suff=filter_suffix)

            # If user requires than unsorted bam be removed:
            if "del_unfiltered" in self.params.keys():
                self.script += "###########\n# Removing unfiltered BAM\n#----------------\n"
                self.script += "\n\nrm -rf %s\n\n" % (active_file)

            # Storing filtered and unfiltered bams:
            self.sample_data["unfiltered_bam"] = self.base_dir + os.path.basename(active_file)  # bam_name
            active_file = active_file + filter_suffix
            self.sample_data["bam"] = self.base_dir + os.path.basename(active_file)
            self.stamp_file(self.sample_data["bam"])

            # The following is so that sort will work on the filtered file without playing around with the sort code:
            bam_name = filtered_name

        if "sort" in self.params.keys():

            self.script += "###########\n# Sorting BAM\n#----------------\n"
            self.script += "%s sort \\\n\t" % self.get_script_env_path()
            if self.params["sort"]:
                self.script += "%s \\\n\t" % self.params["sort"]
            self.script += "-o %s \\\n\t" % (use_dir + os.path.basename(active_file) + sort_suffix)
            self.script += "%s\n\n" % active_file  # (bam_name)

            # If user requires than unsorted bam be removed:
            if "del_unsorted" in self.params.keys():
                self.script += "###########\n# Removing unsorted BAM\n#----------------\n"
                self.script += "\n\nrm -rf %s\n\n" % (active_file)

            # Storing sorted bam in 'bam' slot and unsorted bam in unsorted_bam slot
            self.sample_data["unsorted_bam"] = self.base_dir + os.path.basename(active_file)
            active_file = use_dir + os.path.basename(active_file) + sort_suffix
            self.sample_data["bam"] = self.base_dir + os.path.basename(active_file)
            self.stamp_file(self.sample_data["bam"])

        if "index" in self.params.keys():
            self.script += """
###########
# Indexing BAM
#----------------
{env_script_path} index {params} \\
    {input} \\
    {output}
""".format(env_script_path=self.get_script_env_path(),
                                 params=self.params["index"] if self.params["index"] else "",
                                 input=active_file,
                                 output=use_dir + os.path.basename(active_file) + index_suffix)
            self.sample_data["index"] = self.base_dir + os.path.basename(active_file) + index_suffix
            self.stamp_file(self.sample_data["index"])

        if "flagstat" in self.params.keys():
            self.script += """
###########
#  Calculating BAM statistics
#----------------
{env_script_path} flagstat \\
    {input} \\
    > {output}
""".format(env_script_path=self.get_script_env_path(),
           input=active_file,
           output=use_dir + os.path.basename(active_file) + ".flagstat")

            self.sample_data["flagstat"] = "%s%s.flagstat" % (self.base_dir, os.path.basename(active_file))
            self.stamp_file(self.sample_data["flagstat"])

        if "stats" in self.params.keys():
            self.script += "###########\n# Calculating BAM statistics:\n#----------------\n"
            self.script += "%s stats \\\n\t" % self.get_script_env_path()
            if self.params["stats"]:  # Adding parameters the user might pass
                self.script += "%s \\\n\t" % self.params["stats"]
            self.script += "%s \\\n\t" % active_file  # (use_dir + bam_name)
            self.script += "> %s.stats \n\n" % active_file  # (use_dir + bam_name)
            self.sample_data["stats"] = "%s%s.stats" % (self.base_dir, os.path.basename(active_file))
            self.stamp_file(self.sample_data["stats"])

        if "idxstats" in self.params.keys() and "index" in self.sample_data:
            self.script += "###########\n# Calculating index statistics (idxstats):\n#----------------\n"
            self.script += "%s idxstats \\\n\t" % self.get_script_env_path()
            # idxstats has no uder defined parameters...
            self.script += "%s \\\n\t" % active_file  # (use_dir + bam_name)
            self.script += "> %s.idxstat.tab \n\n" % active_file  # (use_dir + bam_name)
            self.sample_data["idxstats"] = "%s%s.idxstat.tab" % (
                self.base_dir, os.path.basename(active_file))  # bam_name)
            self.stamp_file(self.sample_data["idxstats"])

        # Adding code for fastq or fasta extraction from bam:
        for type in (set(self.params.keys()) & set(["fasta", "fastq"])):
            self.script += """###########\n# Extracting fastq files from BAM:\n#----------------\n"""
            self.script += self.get_script_env_path()
            self.script += "{command} {params} \\\n\t".format(command=type, params=self.params[type])
            if "fastq.F" in self.sample_data:
                self.script += "-1  {readsF} \\\n\t".format(readsF=(active_file + ".F." + type))
                self.script += "-2  {readsR} \\\n\t".format(readsR=(active_file + ".R." + type))
            else:
                self.script += "-0  {readsS} \\\n\t".format(readsS=(active_file + ".S." + type))
            # -s and mixed paired-single not supported yet
            self.script += "{bam} \n\n".format(bam=(active_file))

            # Storing and Stamping files
            if "fastq.F" in self.sample_data:
                self.sample_data["%s.F" % type] = "%s%s.F.%s" % (
                    self.base_dir, os.path.basename(active_file), type)
                self.sample_data["%s.R" % type] = "%s%s.R.%s" % (
                    self.base_dir, os.path.basename(active_file), type)
                self.stamp_file(self.sample_data["%s.F" % type])
                self.stamp_file(self.sample_data["%s.R" % type])
            else:
                self.sample_data["%s.S" % type] = "%s%s.S.%s" % (
                    self.base_dir, os.path.basename(active_file), type)
                self.stamp_file(self.sample_data["%s.S" % type])

        if "del_sam" in self.params.keys() and "sam" in self.sample_data:
            self.script += "###########\n# Removing SAM\n#----------------\n\n"
            self.script += "rm -rf %s\n\n" % self.sample_data["sam"]

        # self.stamp_dir_files(self.base_dir)
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir, self.base_dir)  # Sees to copying local files to final destination (and other stuff)

        self.create_low_level_script()

