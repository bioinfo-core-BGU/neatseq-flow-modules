# -*- coding: UTF-8 -*-
""" 
``RSEM_mapper``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``rsem-calculate-expression``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fasta files in one of the following slots:

    * ``sample_data["project_data"]["fasta.nucl"]``  (``scope`` = ``project``)
    * ``sample_data[<sample>]["fasta.nucl"]``   (``scope`` = ``sample``)
    
* If neither exists, please supply ``reference`` parameter.
       

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Puts output index files in one of the following slot:

    * ``self.sample_data[<sample>]["genes.counts"]``
    * ``self.sample_data[<sample>]["isoforms.counts"]``

And the following BAMs, depending on redirected params:

    * ``self.sample_data[<sample>]["genome.unsorted.bam"]``
    * ``self.sample_data[<sample>]["genome.bam"]``
    * ``self.sample_data[<sample>]["transcript.unsorted.bam"]``
    * ``self.sample_data[<sample>]["transcript.bam"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope", "project | sample", "The scope of the RSEM index. Must match the scope in the RSEM_prep instance."
    "result2use", "genes | isoforms", "Summarize counts at the gene or isoform level."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mapping fastq files::

    RSEM_map:
        module:             RSEM_mapper
        base:               merge1
        script_path:        {Vars.paths.RSEM.rsem-calculate-expression}
        reference:              /path/to/fasta
        redirects:
            --gtf:          /path/to/gtf
            --transcript-to-gene-map: /path/to/map_file

Parsing an existing BAM alignment file::

    RSEM_parse_bam:
        module:         RSEM_mapper
        base:           mv_transcript_bam_to_bam
        script_path:    {Vars.paths.RSEM.rsem-calculate-expression}
        scope:          project
        qsub_params:
            -pe:        shared 20
        redirects:
            --num-threads:  20


References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_RSEM_mapper(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        if "scope" not in self.params:
            raise AssertionExcept("Please supply a scope parameter: either 'sample' or 'project'!")

        if "result2use" not in self.params:
            self.params["result2use"] = "isoforms"
        if self.params["result2use"] not in ["genes","isoforms"]:
            raise AssertionExcept("'result2use' must be either 'genes' or 'isoforms'")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:
                if "RSEM.index" not in self.sample_data[sample]:
                    raise AssertionExcept("No RSEM.index exists for RSEM mapper in sample!", sample)
        else:    #if self.params["scope"] == "project":
            if "RSEM.index" not in self.sample_data["project_data"]:
                raise AssertionExcept("No RSEM.index exists for RSEM mapper!")

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        for sample in self.sample_data["samples"]:

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define location and prefix for output files:
            output_prefix = sample + "_RSEM"
            # Getting alignment file, if exists, in this order: bam, sam, cram
            alignment = None
            for type in ["bam","sam","cram"]:
                if type in self.sample_data[sample]:
                    alignment = self.sample_data[sample][type]
                 
            # Get constant part of script:
            self.script += self.get_script_const()

            if alignment:
                self.script += "--alignments %s \\\n\t" % alignment
            elif "fastq.F" in self.sample_data[sample]:
                self.script += """\
--paired-end \\
\t{readsF} \\
\t{readsR} \\
\t""".format(readsF=self.sample_data[sample]["fastq.F"],
           readsR=self.sample_data[sample]["fastq.R"])
            elif "fastq.S" in self.sample_data[sample]:
                self.script += "{readsS} \\\n\t".format(readsS=self.sample_data[sample]["fastq.S"])
            else:
                raise AssertionExcept("No fastq file. You must have either alignment files or fastq files...\n", sample)

            if self.params["scope"] == "sample":
                self.script += "%s \\\n\t" % self.sample_data[sample]["RSEM.index"]
            else:    # if self.params["scope"] == "project":
                self.script += "%s \\\n\t" % self.sample_data["project_data"]["RSEM.index"]

            self.script += "{dir}{sample}\n\n".format(dir=sample_dir,sample=sample)
            # Saving bam files:
            if "--output-genome-bam" in self.params["redir_params"].keys():
                if "--sort-bam-by-read-name" in self.params["redir_params"] or "--sort-bam-by-coordinate" in self.params["redir_params"]:
                    self.sample_data[sample]["genome.unsorted.bam"] = sample_dir + sample + ".genome.bam"
                    self.sample_data[sample]["genome.bam"] = sample_dir + sample + ".genome.sorted.bam"
                else:
                    self.sample_data[sample]["genome.bam"] = sample_dir + sample + ".genome.bam"
            if "--sort-bam-by-read-name" in self.params["redir_params"] or "--sort-bam-by-coordinate" in self.params["redir_params"]:
                self.sample_data[sample]["transcript.unsorted.bam"] = sample_dir + sample + ".transcript.bam"
                self.sample_data[sample]["transcript.bam"] = sample_dir + sample + ".transcript.sorted.bam"
            else:
                self.sample_data[sample]["transcript.bam"] = sample_dir + sample + ".transcript.bam"
            
            self.sample_data[sample]["genes.counts"] = sample_dir + sample + ".genes.results"
            self.sample_data[sample]["isoforms.counts"] = sample_dir + sample + ".isoforms.results"
            
            self.stamp_file(self.sample_data[sample]["genes.counts"])
            self.stamp_file(self.sample_data[sample]["isoforms.counts"])
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
           
            self.create_low_level_script()

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        results_list = " \\\n\t".join([self.sample_data[sample][self.params["result2use"] + ".counts"]
                                       for sample
                                       in self.sample_data["samples"]])
        output_basename = "{name}.{type}.matrix".format(name=self.sample_data["Title"],
                                                        type=self.params["result2use"])

        self.script = """
{path}{sep}rsem-generate-data-matrix \\
    {result_list} \\
    > {output}

        """.format(path=os.path.dirname(self.params["script_path"]),
                   sep=os.sep,
                   result_list=results_list,
                   output=use_dir + output_basename)
        self.sample_data["project_data"][self.params["result2use"] + ".matrix"] = self.base_dir + output_basename
        self.stamp_file(self.sample_data["project_data"][self.params["result2use"] + ".matrix"])

        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir, self.base_dir)  # Sees to copying local files to final destination (and other stuff)

