# -*- coding: UTF-8 -*-
""" 
``fastq_screen``
-----------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for executing ``fastq_screen`` on sequence files.

Input files are specified with the ``type`` parameter or taken from the fastq slots, one script per fastq file.

In regular mode, no output file are produced. However, if the ``--tag`` is included, the tagged file will be stored in the equivalent ``fastq.X`` slot.
If a ``--filter`` tag is included, the filtered file will be stored in the equivalent ``fastq.X`` slot.

The parameters can be passed through a configuration file specified in the redirected parameters with the ``--conf`` parameter.

Alternatively, if you do not specify the configuration file, one will be produced for you. For this, you must include:

1. A ``genomes`` section specifying genome indices to screen against (see examples below) and
2. an ``aligner`` section specifying the aligning program to use and it's path.

Additionally, if a ``--threads`` parameter is included in the redirects, it will be incorporated into the configuration file.

.. Attention:: If a ``--bisulfite`` redirected parameter is included, it should contain the path to ``Bismark``, which will be included in the configuration file.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files in at least one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* If ``--tag`` and/or ``--filter`` or ``--nohits`` are included, puts output fastq files in:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "genomes", "``name: index`` pairs (see examples) ", "If ``--conf`` not provided, genomes to screen against."
    "aligner", "``name: index`` single pair", "If ``--conf`` not provided, path to aligner to use."
    

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No configuration file::

    fastq_screen:
        module:         fastq_screen
        base:           merge1
        script_path:    {Vars.paths.fastq_screen}
        qsub_params:
            -pe:        shared 60
        aligner:
            bowtie2:    {Vars.paths.bowtie2}
        genomes:
            Human:      {Vars.databases.human}
            Mouse:      {Vars.databases.moiuse}
            PhiX:       {Vars.databases.phix}
        redirects:
            --filter:   200
            --tag:
            # --nohits:
            --force: 
            --threads:  60 

With configuration file::

    fastq_screen:
        module:         fastq_screen
        base:           merge1
        script_path:    {Vars.paths.fastq_screen}
        qsub_params:
            -pe:        shared 60
        redirects:
            --conf:     {Vars.paths.fastq_screen_conf_file}
            --filter:   200
            --tag:
            # --nohits:
            --force: 

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wingett, S.W. and Andrews, S., 2018. **FastQ Screen: A tool for multi-genome mapping and quality control**. F1000Research, 7.

"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_fastq_screen(Step):

    auto_redirs = "--outdir --version --help --get_genomes --aligner".split(" ")
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        if "--conf" not in self.params["redir_params"]:
            if "genomes" not in self.params:
                raise AssertionExcept("Please specify either a configuration file via `--conf` (in redirects) or a "
                                      "list of genomes via `genomes`")
            if "aligner" not in self.params:
                raise AssertionExcept("Please specify either a configuration file via `--conf` (in redirects) or a "
                                      "mapper to use with 'aligner")
            if not isinstance(self.params["aligner"], dict) or len(list(self.params["aligner"].keys())) != 1:
                raise AssertionExcept("'aligner' must be a dictionary containing the path to 'bowtie2', 'bowtie' or 'bwa'")
            if len(list(set(self.params["aligner"].keys()) & set(["bowtie2", "bowtie", "bwa"]))) != 1:
                raise AssertionExcept("'--aligner' must contain one and only one of 'bowtie2', 'bowtie' and 'bwa'")

        if "--filter" in self.params["redir_params"] and "--tag" not in self.params["redir_params"]:
            self.write_warning("When asking for filtered sequences, you must also ask to tag them. Adding --tag.", admonition="attention")
            self.params["redir_params"]["--tag"] = None
        if "--filter" in self.params["redir_params"] or "--tag" in self.params["redir_params"]:
            if "--nohits" in self.params["redir_params"]:
                raise AssertionExcept("Option --nohits may not be specified with --filter/--tag.")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # If --conf not specified, create it:
        if "--conf" not in self.params["redir_params"]:
            self.params["redir_params"]["--conf"] = self.create_conf_file()

        if "scope" not in self.params:
            self.params["scope"] = "sample"
            self.write_warning("No scope defined. Using sample scope by default.", admonition="attention".upper())

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        if "type" in self.params:
            # Convert string to list:
            if isinstance(self.params["type"], str):
                self.params["type"] = [self.params["type"]]
            elif isinstance(self.params["type"],list):
                self.params["type"] = list(set(self.params["type"]))
            else:
                raise AssertionExcept("'type' must be one of fastq.F, fastq.R, fastq.S, or a list thereof")
            # Check legitimacy:
            if any([type not in ["fastq.F", "fastq.R", "fastq.S"] for type in self.params["type"]]):
                raise AssertionExcept("'type' must be one of fastq.F, fastq.R, fastq.S, or a list thereof")
            for sample in sample_list:  # Getting list of samples out of samples_hash
                if any([type not in self.sample_data[sample] for type in self.params["type"]]):
                    raise AssertionExcept(
                        "Sample does not contain file of specified 'type' ({type})".format(type=self.params["type"]),sample)
        else:
            self.params["type"] = ["fastq.F", "fastq.R", "fastq.S"]

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

            for intype in self.params["type"]:
                # Allow different samples to have different types:
                if intype not in self.sample_data[sample]:
                    continue
                infile = self.sample_data[sample][intype]
                outfn = None
                if "--tag" in self.params["redir_params"]:
                    outfn = ".tagged".join(os.path.splitext(os.path.basename(infile)))
                if "--filter" in self.params["redir_params"] or "--nohits" in self.params["redir_params"]:
                    outfn = ".tagged_filter".join(os.path.splitext(os.path.basename(infile)))

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name("_".join([sample,intype]))
                self.script = self.get_script_const()
                self.script += """\
--outdir  {outdir} \\
\t{fastq}
            """.format(fastq=infile, outdir=use_dir)


                self.sample_data[sample][intype+".fastq_screen"] = sample_dir
                if outfn:   # An output was requested
                    self.sample_data[sample][intype] = sample_dir + outfn
                    self.stamp_file(self.sample_data[sample][intype])

                self.local_finish(use_dir, sample_dir)
                self.create_low_level_script()




    def create_conf_file(self):
        """ Build a conf file
        """

        conf_fn = self.base_dir + "fastq_screen.conf"

        with open(conf_fn,"w") as conf:
            conf.write("## Conf file created by NeatSeq-Flow\n\n")
            conf.write("{mapper} {mapper_path}\n\n".format(mapper=list(self.params["aligner"].keys())[0].upper(),
                                                           mapper_path=list(self.params["aligner"].values())[0]))
            if "--bisulfite" in self.params["redir_params"]:
                conf.write("BISMARK {mapper_path}\n\n".format(mapper_path=self.params["redir_params"]["--bisulfite"]))
                self.params["redir_params"]["--bisulfite"] = ""
            if "--threads" in self.params["redir_params"]:
                conf.write("THREADS {threads}\n\n".format(threads=self.params["redir_params"]["--threads"]))
            for genome in self.params["genomes"]:
                conf.write("DATABASE	{name}	{path}\n\n".format(name=genome,path=self.params["genomes"][genome]))

        conf.close()
        return conf_fn

